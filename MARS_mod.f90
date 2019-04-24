! <MARS_mod.f90 - A component of the EMEP MSC-W Unified Eulerian
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
 
module MARS_mod
 ! -----------------------------------------------------------------------
 ! Calculates gas-aerosol equilibrium for SO4, HNO3-NO3 and NH3-NH4 system
 ! Made available by Frank Binkowski (originally from EPA's RPM model)
 ! Presently not in use, EQSAM is used for inorganic equilibrium stuff
 !------------------------------------------------------------------------

 use CheckStop_mod,       only: CheckStop
 use Io_mod,              only: ios, datewrite
 use MARS_Aero_water_mod, only:  Awater!,awater_2900
 use Config_module,       only: NPROC
 use Debug_module,        only: DEBUG   ! -> DEBUG%EQUIB
 use Par_mod,             only: me      ! processor number
 implicit none
 private

!to test without the smoothing between different TNH4/TSO4 regimes
!before rv4.2beta (Jan 2013), the code was (or should be) equivalent to MARS_RATIO_SMOOTH=.false.
  logical, parameter :: MARS_RATIO_SMOOTH=.true.


  real, parameter ::    FLOOR = 1.0E-30       ! minimum concentration  
                                              ! -30 from RPM


 !/- subroutines:
 public   ::  rpmares,   &
              rpmares_2900,   &!svn version 2908, for testing if there are significant differences. will be deleted
              rpmares_new,   &!GEOSCHEM adapted version
              DO_RPMARES_new, &!driver for rpmares_new
              awater_new,   &!almost identical to awater. Only difference Y2     = ( 1.0d0 - MFSSO4 ) / MFSSO4?
              cubic,     &
              actcof,    &
              IS_SAFE_DIV

      integer, private, save :: MAXNNN1 = 0
      integer, private, save :: MAXNNN2 = 0
!ds      real        MWNO3            ! molecular weight for NO3
      real, private, parameter :: MWNO3  = 62.0049 

!ds      real        MWHNO3           ! molecular weight for HNO3
      real, private, parameter :: MWHNO3 = 63.01287       

!ds      real        MWSO4            ! molecular weight for SO4
      real, private, parameter :: MWSO4 = 96.0576 

!ds      real        MWHSO4           ! molecular weight for HSO4
      real, private, parameter :: MWHSO4 = MWSO4 + 1.0080

!ds      real        MH2SO4           ! molecular weight for H2SO4
      real, private, parameter :: MH2SO4 = 98.07354 

!ds      real        MWNH3            ! molecular weight for NH3
      real, private, parameter :: MWNH3 = 17.03061 

!ds      real        MWNH4            ! molecular weight for NH4
      real, private, parameter :: MWNH4 = 18.03858
 contains

      real function poly4 (a,x)

! Calculates the polynomial based on 4 coefficients a(1:4):

!-- arguments
      real, dimension(4), intent(in)  ::  a
      real, intent(in)  ::  x

      poly4 = a(1) + x * ( a(2) + x * ( a(3) + x * ( a(4) ) ) )
      
      end function poly4

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      real function poly6(a,x)

! Calculates the polynomial based on 6 coefficients a(1:6):

!-- arguments
      real, dimension(6), intent(in)  ::  a
      real, intent(in)  ::  x

      poly6 = a(1) + x * ( a(2) + x * ( a(3) + x * ( a(4) +  &
              x * ( a(5) + x * (a(6)  ) ) ) ) )

      end  function poly6 


!from GEOS-Chem
!Could/should be moved to another module
     FUNCTION IS_SAFE_DIV( N, D, R4 ) RESULT( F )
!
! !INPUT PARAMETERS: 
!
      REAL*8,  INTENT(IN)           :: N    ! Numerator
      REAL*8,  INTENT(IN)           :: D    ! Denominator
      LOGICAL, INTENT(IN), OPTIONAL :: R4   ! Logical flag to use the limits 
                                            !  of REAL*4 to define underflow
                                            !  or overflow.  Extra defensive.
!
! !OUTPUT PARAMETERS:
!
      LOGICAL                       :: F    ! =F if division isn't allowed
                                            ! =T otherwise
!
! !REMARKS:
!  UnderFlow, OverFlow and NaN are tested for. If you need to
!  differentiate between the three, use the SAFE_DIV (phs, 4/14/09)
!
! !REVISION HISTORY:
!  11 Jun 2008 - P. Le Sager - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      INTEGER MaxExp, MinExp
      REAL*4  RR

      !==================================================================
      ! IS_SAFE_DIV begins here!
      !==================================================================

      MaxExp = MAXEXPONENT( N )
      MinExp = MINEXPONENT( N )

      IF ( PRESENT( R4 ) ) THEN
         IF ( R4 ) THEN
            MaxExp = MAXEXPONENT( RR )
            MinExp = MINEXPONENT( RR )
         ENDIF
      ENDIF  

      IF ( EXPONENT(N) - EXPONENT(D) >= MaxExp .or. D==0 .or.&
          EXPONENT(N) - EXPONENT(D) <= MinExp  ) THEN
         F = .FALSE.
      ELSE
         F = .TRUE.
      ENDIF

      ! Return to calling program
      END FUNCTION IS_SAFE_DIV

!>-------------------------------------------------------------------------------<
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      subroutine rpmares ( SO4, HNO3, NO3, NH3, NH4, RH, TEMP,   &
                           ASO4, ANO3, AH2O, ANH4, GNH3, GNO3,   &
                           ERRMARK,debug_flag) 

!-----------------------------------------------------------------------
!C
!C Description:
!C
!C   ARES calculates the chemical composition of a sulfate/nitrate/
!C   ammonium/water aerosol based on equilibrium thermodynamics.
!C
!C   This code considers two regimes depending upon the molar ratio 
!C   of ammonium to sulfate. 
!C
!C   For values of this ratio less than 2,the code solves a cubic for 
!C   hydrogen ion molality, HPLUS,  and if enough ammonium and liquid
!C   water are present calculates the dissolved nitric acid. For molal
!C   ionic strengths greater than 50, nitrate is assumed not to be present. 
!C   
!C   For values of the molar ratio of 2 or greater, all sulfate is assumed
!C   to be ammonium sulfate and a calculation is made for the presence of
!C   ammonium nitrate.
!C
!C   The Pitzer multicomponent approach is used in subroutine ACTCOF to
!C   obtain the activity coefficients. Abandoned -7/30/97 FSB 

!c   The Bromley method of calculating the activity coefficients is s used
!c    in this version

!c   The calculation of liquid water
!C   is done in subroutine water. Details for both calculations are given
!C   in the respective subroutines.
!C
!C   Based upon MARS due to 
!C   P. Saxena, A.B. Hudischewskyj, C. Seigneur, and J.H. Seinfeld, 
!C   Atmos. Environ., vol. 20, Number 7, Pages 1471-1483, 1986.
!C
!C   and SCAPE due to 
!C   Kim, Seinfeld, and Saxeena, Aerosol Ceience and Technology,
!C   Vol 19, number 2, pages 157-181 and pages 182-198, 1993.
!C
!C NOTE: All concentrations supplied to this subroutine are TOTAL
!C       over gas and aerosol phases
!C
!C Parameters:
!C 
!C  SO4   : Total sulfate in MICROGRAMS/M**3 as sulfate (IN)
!C  HNO3  : Nitric Acid in MICROGRAMS/M**3 as nitric acid (IN)
!C  NO3   : Total nitrate in MICROGRAMS/M**3 as nitric acid (IN)
!C  NH3   : Total ammonia in MICROGRAMS/M**3 as ammonia (IN)
!C  NH4   : Ammonium in MICROGRAMS/M**3 as ammonium (IN)
!C  RH    : Fractional relative humidity (IN)
!C  TEMP  : Temperature in Kelvin (IN)
!C  GNO3  : Gas phase nitric acid in MICROGRAMS/M**3 (OUT)
!C  GNH3  : Gas phase ammonia in MICROGRAMS/M**3 (OUT)
!C  ASO4  : Aerosol phase sulfate in MICROGRAMS/M**3 (OUT) 
!C  ANO3  : Aerosol phase nitrate in MICROGRAMS/M**3 (OUT)
!C  ANH4  : Aerosol phase ammonium in MICROGRAMS/M**3 (OUT)
!C  AH2O  : Aerosol phase water in MICROGRAMS/M**3 (OUT)
!C  NITR  : Number of iterations for obtaining activity coefficients  (OUT) 
!C  NR    : Number of real roots to the cubic in the low ammonia case (OUT)
!C 
!C Revision History:
!C      Who       When        Detailed description of changes
!C   ---------   --------  -------------------------------------------
!C   S.Roselle   11/10/87  Received the first version of the MARS code
!C   S.Roselle   12/30/87  Restructured code
!C   S.Roselle   2/12/88   Made correction to compute liquid-phase 
!C                         concentration of H2O2.  
!C   S.Roselle   5/26/88   Made correction as advised by SAI, for 
!C                         computing H+ concentration.
!C   S.Roselle   3/1/89    Modified to operate with EM2
!C   S.Roselle   5/19/89   Changed the maximum ionic strength from 
!C                         100 to 20, for numerical stability.
!C   F.Binkowski 3/3/91    Incorporate new method for ammonia rich case
!C                         using equations for nitrate budget.
!C   F.Binkowski 6/18/91   New ammonia poor case which 
!C                         omits letovicite.
!C   F.Binkowski 7/25/91   Rearranged entire code, restructured
!C                         ammonia poor case.
!C   F.Binkowski 9/9/91    Reconciled all cases of ASO4 to be output
!C                         as SO4--
!C   F.Binkowski 12/6/91   Changed the ammonia defficient case so that 
!C                         there is only neutralized sulfate (ammonium
!C                         sulfate) and sulfuric acid.
!C   F.Binkowski 3/5/92    Set RH bound on AWAS to 37 % to be in agreement 
!C                          with the Cohen et al. (1987)  maximum molality
!C                          of 36.2 in Table III.( J. Phys Chem (91) page
!C                          4569, and Table IV p 4587.)
!C   F.Binkowski 3/9/92    Redid logic for ammonia defficient case to remove
!C                         possibility for denomenator becoming zero; 
!C                         this involved solving for HPLUS first.
!C                         Note that for a relative humidity
!C                          less than 50%, the model assumes that there is no 
!C                          aerosol nitrate.
!C   F.Binkowski 4/17/95   Code renamed  ARES (AeRosol Equilibrium System)  
!C                          Redid logic as follows
!C                         1. Water algorithm now follows Spann & Richardson
!C                         2. Pitzer Multicomponent method used
!C                         3. Multicomponent practical osmotic coefficient 
!C                            use to close iterations.
!C                         4. The model now assumes that for a water
!C                            mass fraction WFRAC less than 50% there is
!C                            no aerosol nitrate.
!C   F.Binkowski 7/20/95   Changed how nitrate is calculated in ammonia poor 
!C                         case, and changed the WFRAC criterion to 40%.
!C                         For ammonium to sulfate ratio less than 1.0 
!C                         all ammonium is aerosol and no nitrate aerosol
!C                         exists.
!C   F.Binkowski 7/21/95   Changed ammonia-ammonium in ammonia poor case to
!C                         allow gas-phase ammonia to exist. 
!C   F.Binkowski 7/26/95   Changed equilibrium constants to values from 
!C                         Kim et al. (1993) 
!C   F.Binkowski 6/27/96   Changed to new water format
!c   F.Binkowski 7/30/97   Changed to Bromley method for multicomponent 
!c                         activity coefficients. The binary activity coefficients 
!c                         are the same as the previous version
!c   F.Binkowski 8/1/97    Chenged minimum sulfate from 0.0 to 1.0e-6 i.e.
!c                         1 picogram per cubic meter
!C   I.Ackermann 2/23/98   modification for solving precision problems
!C        	  on workstations
!C   I.Ackermann 2/25/99   changed logic as follows:
!c                         If iterations fail, initial values of nitrate
!c                          are retained.
!c                         Total mass budgets are changed to account for gas
!c                         phase returned. (incorporated from FZB's models3
!c                         framework)
!C                         eliminated ratio=5 !! for low to zero sulfate
!C   I.Ackermann 3/17/99   modified ratio = 5 treatment see RB3,p.125
!C
!C-----------------------------------------------------------------------


!...........ARGUMENTS and their descriptions

  real, intent(in) ::  SO4   &     ! Total sulfate in micrograms / m**3 
                      ,HNO3  &     ! Total nitric acid in micrograms / m**3
                      ,NO3   &     ! Total nitrate in micrograms / m**3
                      ,NH3   &     ! Total ammonia in micrograms / m**3
                      ,NH4   &     ! Total ammonium in micrograms / m**3
                      ,RH    &     ! Fractional relative humidity 
                      ,TEMP        ! Temperature in Kelvin 

  real, intent(out)::  ASO4  &     ! Aerosol sulfate in micrograms / m**3 
                      ,ANO3  &     ! Aerosol nitrate in micrograms / m**3
                      ,AH2O  &     ! Aerosol liquid water content water in micrograms / m**3
                      ,ANH4  &     ! Aerosol ammonium in micrograms / m**3
                      ,GNO3  &     ! Gas-phase nitric acid in micrograms / m**3
                      ,GNH3        ! Gas-phase ammonia in micrograms / m**3

  logical, intent(in) :: debug_flag

!C...........INCLUDES and their descriptions
!!      INCLUDE SUBST_CONST          ! constants

!...........PARAMETERS and their descriptions:

      real        MWNACL           ! molecular weight for NaCl
      parameter ( MWNACL = 58.44277 )

!ds moved a bunch upstairs.

      real        MWORG            ! molecular weight for Organic Species
      parameter ( MWORG = 16.0 )

      real        MWCL             ! molecular weight for Chloride  
      parameter ( MWCL = 35.453 )

      real        MWAIR            ! molecular weight for AIR
      parameter ( MWAIR = 28.964 )

      real        MWLCT            ! molecular weight for Letovicite
      parameter ( MWLCT = 3.0 * MWNH4 + 2.0 * MWSO4 + 1.0080 )

      real        MWAS             ! molecular weight for Ammonium Sulfate
      parameter ( MWAS = 2.0 * MWNH4 + MWSO4 )

      real        MWABS            ! molecular weight for Ammonium Bisulfate 
      parameter ( MWABS = MWNH4 + MWSO4 + 1.0080 )


!...........SCRATCH LOCAL VARIABLES and their descriptions:
       
      REAL        fRH              ! Index set to percent relative humidity  
      INTEGER     NITR             ! Number of iterations for activity coefficients
      INTEGER     NNN              ! Loop index for iterations 
      INTEGER     NR               ! Number of roots to cubic equation for HPLUS
      INTEGER     ERRMARK

      real          A0             ! Coefficients and roots of 
      real          A1             ! Coefficients and roots of 
      real          A2             ! Coefficients and roots of 
      real        AA               ! Coefficients and discriminant for quadratic equation for ammonium nitrate
      real        BAL              ! internal variables ( high ammonia case)
      real        BB               ! Coefficients and discriminant for quadratic equation for ammonium nitrate
      real        BHAT             ! Variables used for ammonia solubility 
      real        CC               ! Coefficients and discriminant for quadratic equation for ammonium nitrate
      real        CONVT            ! Factor for conversion of units  
      real        DD               ! Coefficients and discriminant for quadratic equation for ammonium nitrate
      real        DISC             ! Coefficients and discriminant for quadratic equation for ammonium nitrate
      real        EROR             ! Relative error used for convergence test  
      real        FNH3             ! "Free ammonia concentration", that which exceeds TWOSO4       
      real        GAMAAB           ! Activity Coefficient for (NH4+, HSO4-)GAMS( 2,3 )
      real        GAMAAN           ! Activity coefficient for (NH4+, NO3-) GAMS( 2,2 )
      real        GAMAHAT          ! Variables used for ammonia solubility 
      real        GAMANA           ! Activity coefficient for (H+ ,NO3-)   GAMS( 1,2 )
      real        GAMAS1           ! Activity coefficient for (2H+, SO4--) GAMS( 1,1 )
      real        GAMAS2           ! Activity coefficient for (H+, HSO4-)  GAMS( 1,3 )
      real        GAMOLD           ! used for convergence of iteration
      real        GASQD            ! internal variables ( high ammonia case)
      real        HPLUS            ! Hydrogen ion (low ammonia case) (moles / kg water)
      real        K1A              ! Equilibrium constant for ammoniua to ammonium
      real        K2SA             ! Equilibrium constant for sulfate-bisulfate (aqueous)
      real        K3               ! Dissociation constant for ammonium nitrate 
      real        KAN              ! Equilibrium constant for ammonium nitrate (aqueous)
      real        KHAT             ! Variables used for ammonia solubility 
      real        KNA              ! Equilibrium constant for nitric acid (aqueous)   
      real        KPH              ! Henry's Law Constant for ammonia       
      real        KW               ! Equilibrium constant for water dissociation             
      real        KW2              ! Internal variable using KAN 
      real        MAN              ! Nitrate (high ammonia case) (moles / kg water)
      real        MAS              ! Sulfate (high ammonia case) (moles / kg water)
      real        MHSO4            ! Bisulfate (low ammonia case) (moles / kg water)
      real        MNA              ! Nitrate (low ammonia case) (moles / kg water)
      real        MNH4             ! Ammonium (moles / kg water)
      real        MOLNU            ! Total number of moles of all ions
      real        MSO4             ! Sulfate (low ammonia case) (moles / kg water)
      real        PHIBAR           ! Practical osmotic coefficient      
      real        PHIOLD           ! Previous value of practical osmotic coefficient used for convergence of iteration
      real        RATIO            ! Molar ratio of ammonium to sulfate
      real        RK2SA            ! Internal variable using K2SA
      real        RKNA             ! Internal variables using KNA
      real        RKNWET           ! Internal variables using KNA
      real        RR1
      real        RR2
      real        STION            ! Ionic strength
      real        T1               ! Internal variables for temperature corrections
      real        T2               ! Internal variables for temperature corrections
      real        T21              ! Internal variables of convenience (low ammonia case)
      real        T221             ! Internal variables of convenience (low ammonia case)
      real        T3               ! Internal variables for temperature corrections
      real        T4               ! Internal variables for temperature corrections
      real        T6               ! Internal variables for temperature corrections
      real        TNH4             ! Total ammonia and ammonium in micromoles / meter ** 3
      real        TNO3             ! Total nitrate in micromoles / meter ** 3
      real        TOLER1           ! Tolerances for convergence test 
      real        TOLER2           ! Tolerances for convergence test 
      real        TSO4             ! Total sulfate in micromoles / meter ** 3
      real        TWOSO4           ! 2.0 * TSO4  (high ammonia case) (moles / kg water)
      real        WFRAC            ! Water mass fraction 
      real        WH2O             ! Aerosol liquid water content (internally)
                                   ! micrograms / meter **3 on output
                                   ! internally it is 10 ** (-6) kg (water) / meter ** 3
                                   ! the conversion factor (1000 g = 1 kg) is applied 
                                   ! for AH2O output
      real        WSQD             ! internal variables ( high ammonia case)
      real        XNO3             ! Nitrate aerosol concentration in micromoles / meter ** 3
      real        XXQ              ! Variable used in quadratic solution
      real        YNH4             ! Ammonium aerosol concentration in micromoles / meter** 3
      real        ZH2O             ! Water variable saved in case ionic strength too high.
      real        ZSO4             ! Total sulfate molality - mso4 + mhso4 (low ammonia case) (moles / kg water)

      real        CAT( 2 )         ! Array for cations (1, H+); (2, NH4+) (moles / kg water)
      real        AN ( 3 )         ! Array for anions (1, SO4--); (2, NO3-); (3, HSO4-)  (moles / kg water) 
      real        CRUTES( 3 )      ! Coefficients and roots of 
      real        GAMS( 2, 3 )     ! Array of activity coefficients 
      real        MINSO4           ! Minimum value of sulfate laerosol concentration
       parameter( MINSO4 = 1.0E-6 / MWSO4 ) 
      real        MINNO3
       parameter( MINNO3 = 1.0E-6 / MWNO3 )    !2/25/99 IJA
!st      real        FLOOR
!st       parameter( FLOOR = 1.0E-30) ! minimum concentration       
!2/25/99 IJA
! FSB New variables Total ammonia and nitrate mass concentrations
      real  TMASSNH3  ! Total ammonia (gas and particle)
                      !  as ammonia gas [ug /m**3]
      real  TMASSHNO3 ! Total nitrate (vapor and particle) as
                      !  nitric acid vapor [ug /m**3]

!for RATIO between RATIO_Low and RATIO_High, consider the species in two boxes one with
! RATIO_Low and the other with  RATIO_High. Then take a linear combination of the results in each box
      logical ::DO_RATIO_High_2,DO_RATIO_Low_2
      real, parameter ::RATIO_Low = 1.95,RATIO_High= 2.05 !somewhat arbitrarily
      real ::ASO4_High, ANO3_High, AH2O_High, ANH4_High, GNH3_High, GNO3_High
      real ::ASO4_High1, ANO3_High1, AH2O_High1, ANH4_High1, GNH3_High1, GNO3_High1,diff1
      real ::ASO4_Low, ANO3_Low, AH2O_Low, ANH4_Low, GNH3_Low, GNO3_Low
      real ::ASO4_Low1, ANO3_Low1, AH2O_Low1, ANH4_Low1, GNH3_Low1, GNO3_Low1,diff2
      real::TSO4_HighA,TSO4_LowA,High_Factor,X
!-----------------------------------------------------------------------
      integer, save :: nMarsErrors = 0 ! tracks solution failures DSMARS
!-----------------------------------------------------------------------
!  begin body of subroutine RPMARES
                                                                         
!ASO4=FLOOR;ANO3=FLOOR;AH2O=FLOOR;ANH4=FLOOR;GNO3=FLOOR;GNH3=FLOOR 
!Initialise the output variables

      ASO4=0.0;ANO3=0.0;AH2O=0.0;ANH4=0.0;GNO3=0.0;GNH3=0.0 

      ASO4 = SO4   !ds from RPM

!...convert into micromoles/m**3
 
!..iamodels3 merge NH3/NH4 , HNO3,NO3 here
      TSO4 = MAX( 0.0, SO4 / MWSO4  )
      TSO4 = MAX( FLOOR, TSO4  )      !GEOS added
      TNO3 = MAX( 0.0, (NO3 / MWNO3 + HNO3 / MWHNO3) )
      TNH4 = MAX( 0.0, (NH3 / MWNH3 + NH4 / MWNH4)  )

!2/25/99 IJA
!      TMASSNH3  = MAX(0.0, NH3 + (MWNH3 / MWNH4)  * NH4 )
!      TMASSHNO3 = MAX(0.0, NO3 + (MWHNO3 / MWNO3) * NO3 )

      TMASSNH3  = MAX(0.0, NH3 +  NH4 )
      TMASSHNO3 = MAX(0.0, HNO3 + NO3 )
 
!...now set humidity index fRH as a percent

!st      IRH = NINT( 100.0 * RH )
         fRH = RH 

!GEOS added
      ! For extremely low relative humidity ( less than 1% ) set the 
      ! water content to a minimum and skip the calculation.
      IF ( fRH .LT. 0.01 ) THEN
         AH2O = FLOOR

!not from GEOS:
          ASO4 = SO4
          ANO3 = NO3   
          WH2O = 0.0
          GNH3 = NH3
          ANH4 = NH4
          GNO3 = HNO3

          RETURN
      ENDIF 


!...Check for valid fRH

       fRH = MAX( 0.01, fRH )
       fRH = MIN( 0.99, fRH )

!...Specify the equilibrium constants at  correct
!...  temperature.  Also change units from ATM to MICROMOLE/M**3 (for KAN,
!...  KPH, and K3 )
!...  Values from Kim et al. (1993) except as noted.
 
      CONVT = 1.0 / ( 0.082 * TEMP ) 
      T6 = 0.082E-9 *  TEMP
      T1 = 298.0 / TEMP
      T2 = ALOG( T1 )
      T3 = T1 - 1.0
      T4 = 1.0 + T2 - T1
      KNA  = 2.511E+06 *  EXP(  29.17 * T3 + 16.83 * T4 ) * T6
      K1A  = 1.805E-05 *  EXP(  -1.50 * T3 + 26.92 * T4 )
      K2SA = 1.015E-02 *  EXP(   8.85 * T3 + 25.14 * T4 )
      KW   = 1.010E-14 *  EXP( -22.52 * T3 + 26.92 * T4 )
      KPH  = 57.639    *  EXP(  13.79 * T3 - 5.39  * T4 ) * T6
!!!      K3   =  5.746E-17 * EXP( -74.38 * T3 + 6.12  * T4 ) * T6 * T6  
      KHAT =  KPH * K1A / KW  
      KAN  =  KNA * KHAT  

!...Compute temperature dependent equilibrium constant for NH4NO3
!...  ( from Mozurkewich, 1993)

      K3 = EXP( 118.87  - 24084.0 / TEMP -  6.025  * ALOG( TEMP ) )

!...Convert to (micromoles/m**3) **2

      K3 = K3 * CONVT * CONVT

      WH2O   = 0.0
      STION  = 0.0
      AH2O   = 0.0
      MAS    = 0.0
      MAN    = 0.0
      HPLUS  = 0.0
      TOLER1 = 0.00001
      TOLER2 = 0.001
      NITR   = 0
      NR     = 0
      RATIO  = 0.0
      GAMAAN = 1.0
      GAMOLD = 1.0

!ds from RPM, but removed FLOOR : (slighty different logic)
      if( (TSO4 <= MINSO4 ) .and. (TNO3 < MINNO3) ) then
        !print *, "DSX  HERE"
          ASO4 = SO4   ! MAX..
          ANO3 = NO3   ! MAX..      
          WH2O = 0.0
          AH2O = 0.0
          GNH3 = NH3   ! MAX(FLOOR,NH3)
          GNO3 = HNO3   ! MAX(FLOOR,NO3)
          ANH4 = NH4
          RETURN
       END IF
!ds end rpm


!...set the ratio according to the amount of sulfate and nitrate

      IF ( TSO4 > MINSO4 ) THEN
        RATIO = TNH4 / TSO4

!...If there is no sulfate and no nitrate, there can be no ammonium
!...  under the current paradigm. Organics are ignored in this version.

      ELSE 
      
        IF ( TNO3 <= MINNO3 ) THEN

! *** If there is very little sulfate and nitrate set concentrations
!      to a very small value and return.    
! Jun 2012, Note these values are set in the initialisation
          ASO4 = SO4 ! MAX(FLOOR, ASO4)
          ANO3 = NO3 ! MAX(FLOOR, ANO3 )          
          WH2O = 0.0
          AH2O = 0.0
          GNH3 = NH3 ! MAX(FLOOR,GNH3)
          GNO3 = HNO3 ! MAX(FLOOR,GNO3)
          ANH4 = NH4
          RETURN
       END IF
       
!...For the case of no sulfate and nonzero nitrate, set ratio to 5
!...  to send the code to the high ammonia case if there is more
!...  ammonia than sulfate, otherwise send to low ammonia case.

       IF (TNH4 > TSO4) THEN
         RATIO = 5.0        !this is a high ammonia case with low sulfur
       ELSE
         RATIO = 1.        !this is a low ammonia case with low sulfur
       END IF
 
      END IF 


      DO_RATIO_High_2=.false.
      DO_RATIO_Low_2=.false.
      ASO4_High = 0.0
      ANO3_High = 0.0
      AH2O_High = 0.0
      ANH4_High = 0.0
      GNH3_High = 0.0
      GNO3_High = 0.0
      ASO4_Low = 0.0
      ANO3_Low = 0.0
      AH2O_Low = 0.0
      ANH4_Low = 0.0
      GNH3_Low = 0.0
      GNO3_Low = 0.0
      TSO4_HighA = TSO4
      TSO4_LowA  = TSO4
       IF ( RATIO > 2.0 )then
         DO_RATIO_High_2=.true.
         High_Factor=1.0
      else
         DO_RATIO_Low_2=.true.
         High_Factor=0.0
      end if

      IF ( RATIO >RATIO_Low  .and.  RATIO < RATIO_High .and. MARS_RATIO_SMOOTH)then
         DO_RATIO_High_2=.true.
         DO_RATIO_Low_2=.true.
         TSO4_HighA=TSO4*Ratio/RATIO_High
         TSO4_LowA=TSO4*Ratio/RATIO_Low
!         High_Factor=(RATIO-RATIO_Low)/(RATIO_High-RATIO_Low)
         High_Factor=(RATIO*RATIO_High-RATIO_Low*RATIO_High)/(RATIO*RATIO_High-RATIO*RATIO_Low)
      end if

!....................................
!......... High Ammonia Case ........
!....................................
 
!      IF ( RATIO > 2.0 ) THEN    ! NH4/SO4 > 2
      IF ( DO_RATIO_High_2 ) THEN    ! NH4/SO4 > 2
        GAMAAN = 0.1

!...Set up twice the sulfate for future use.

        TWOSO4 = 2.0 * TSO4_HighA
        XNO3 = 0.0            
        YNH4 = TWOSO4

!...Treat different regimes of relative humidity 

!...ZSR relationship is used to set water levels. Units are
!...  10**(-6) kg water/ (cubic meter of air)
!...  start with ammomium sulfate solution without nitrate

      CALL awater(fRH,TSO4_HighA,YNH4,TNO3,AH2O ) !**** note TNO3
        WH2O = 1.0E-3 * AH2O  
        ASO4_High = TSO4_HighA * MWSO4
        ANO3_High = 0.0
        ANH4_High = YNH4 * MWNH4
if(debug_flag) call datewrite("MARS debug ", -1,(/ ASO4_High, ANH4_High, AH2O /) )

if ( DEBUG%EQUIB ) then
  if( ASO4_High + ANH4_High +  AH2O < 1.0-10 ) then
     call datewrite("MARS failing? ", -1,(/ ASO4_High, ANH4_High, AH2O /) )
     print *, "MARS PROB ", ASO4_High, ANH4_High, AH2O, TSO4_HighA, YNH4
     call CheckStop("MARS")
  end if
end if
        WFRAC = (AH2O + FLOOR)  / ( ASO4_High + ANH4_High +  AH2O + FLOOR  )
        !ds WFRAC = AH2O / ( ASO4 + ANH4 +  AH2O + FLOOR  )
        !CRUDE FIX? WFRAC = AH2O / ( ASO4 + ANH4 +  AH2O )
!!!!       IF ( WFRAC == 0.0 )  RETURN   ! No water       
        IF ( WFRAC < 0.2 ) THEN
 
!..."dry" ammonium sulfate and ammonium nitrate
!...  compute free ammonia

          FNH3 = TNH4 - TWOSO4
          CC = TNO3 * FNH3 - K3

!...check for not enough to support aerosol      

          !dsjIF ( CC <= 0.0 ) THEN
          IF ( CC < FLOOR ) THEN
            XNO3 = 0.0
          ELSE
            AA = 1.0
            BB = -( TNO3 + FNH3 ) 
            DISC = BB * BB - 4.0 * CC

!...Check for complex roots of the quadratic
!...  set nitrate to zero and RETURN if complex roots are found
!2/25/99 IJA

            !DS IF ( DISC < 0.0 ) THEN
            !dsDSdIF ( DISC < FLOOR ) THEN
            IF ( DISC < 0.0 ) THEN
if( DEBUG%EQUIB .and. debug_flag ) print *, "MARS DISC NEG ", XNO3, WH2O, DISC
              XNO3 = 0.0
              AH2O_High = 1000.0 * WH2O
              YNH4 = TWOSO4
              GNO3_High = HNO3
              ASO4_High = TSO4_HighA * MWSO4
              ANO3_High = NO3
              ANH4_High = YNH4 * MWNH4
              GNH3_High = TMASSNH3 - ANH4
              if( GNH3 < 0.0 ) then
                 print *, " NEG GNH3", TWOSO4, ANH4, TMASSNH3
                 call CheckStop("NEG GNH3")
              end if
!              RETURN
          goto 333
            END IF

!...to get here, BB .lt. 0.0, CC .gt. 0.0 always      

            DD = SQRT( DISC )
            XXQ = -0.5 * ( BB + SIGN ( 1.0, BB ) * DD )

!...Since both roots are positive, select smaller root.      

            XNO3 = MIN( XXQ / AA, CC / XXQ )
          
          END IF
!2/25/99 IJA
          AH2O_High = 1000.0 * WH2O
          YNH4 = TWOSO4 + XNO3
          ASO4_High = TSO4_HighA * MWSO4
          !dsSAFE ANO3 = XNO3 * MWNO3
          ANO3_High = min(XNO3 * MWNO3, TMASSHNO3 )
          !dsSAFE ANH4 = YNH4 * MWNH4 ! ds should be safe as NH4/SO4 >2, but anyway:
          ANH4_High = min(YNH4 * MWNH4, TMASSNH3 )  ! ds should be safe as NH4/SO4 >2, but anyway:
          GNH3_High = TMASSNH3 - ANH4_High
          GNO3_High = TMASSHNO3 - ANO3_High
          !    if( GNH3 < 0.0 .or. GNO3 < 0.0 ) then
          !       print *, " NEG GNH3 GNO3", TWOSO4, ANH4, TMASSNH3, ANO3, TMASSHNO3
          !       call CheckStop("NEG GNH3 GNO3")
          !    end if
!          RETURN
          goto 333
        END IF

!...liquid phase containing completely neutralized sulfate and
!...  some nitrate.  Solve for composition and quantity.
 
        MAS = TSO4_HighA / max(WH2O,FLOOR)
        MAN = 0.0
        XNO3 = 0.0
        YNH4 = TWOSO4
        PHIOLD = 1.0

!...Start loop for iteration
 
!...The assumption here is that all sulfate is ammonium sulfate,
!...  and is supersaturated at lower relative humidities.
        diff1=1.0/TOLER1

        DO 1501 NNN = 1, 150 
          NITR = NNN
          GASQD = GAMAAN * GAMAAN
          WSQD = WH2O * WH2O
          KW2 = KAN * WSQD / GASQD
          AA = 1.0 - KW2
          BB = TWOSO4 + KW2 * ( TNO3 + TNH4 - TWOSO4 )
          CC = -KW2 * TNO3 * ( TNH4 - TWOSO4 )

!...This is a quadratic for XNO3 [MICROMOLES / M**3] of nitrate in solution.

          DISC = BB * BB - 4.0 * AA * CC

          if( DEBUG%EQUIB ) then
            MAXNNN1 = NNN
            !if( MAXNNN1 > 140) print "(a,i4,9es12.3)", "NNN1 ", NNN, DISC, TNO3, TNH4, TWOSO4
          end if
!...Check for complex roots, retain inital values and RETURN
!2/25/99 IJA

          !DS IF ( DISC < 0.0 ) THEN
          !dsDS IF ( DISC < FLOOR ) THEN
          IF ( DISC < 0.0 ) THEN
if( DEBUG%EQUIB .and. debug_flag ) print *, "MARS DISC NEG2 ", XNO3, WH2O, DISC
            XNO3 = 0.0
            AH2O_High = 1000.0 * WH2O
            YNH4 = TWOSO4
            GNO3_High = HNO3
            ASO4_High = TSO4_HighA * MWSO4
            ANO3_High = NO3
            !ANH4 = YNH4 * MWNH4
            ANH4_High = min( YNH4 * MWNH4, TMASSNH3)  ! ds added "min"
            GNH3_High = TMASSNH3 - ANH4_High

              !WRITE( 10, * ) ' COMPLEX ROOTS '
!            RETURN
                goto 333
          END IF
! 2/25/99 IJA

! Deal with degenerate case (yoj)

          !DS IF ( AA /= 0.0 ) THEN
          IF ( abs(AA) > FLOOR  ) THEN
             if( DEBUG%EQUIB .and. debug_flag ) print "(a,9es11.3)", "MARS DEGEN  ",  XNO3, WH2O, DISC, AA, BB, CC
             DD = SQRT( DISC )
             XXQ = -0.5 * ( BB + SIGN ( 1.0, BB ) * DD )
             RR1 = XXQ / AA

             RR2 = CC / XXQ! Normally: BB>0, xxQ<0

             
             !...choose minimum positve root         
             
             IF ( ( RR1 * RR2 ) < 0.0 ) THEN
                if( DEBUG%EQUIB .and. debug_flag ) print "(a,10es10.3)", "MARS RR1*RR2  ", XNO3, WH2O, DISC, RR1, RR2
                XNO3 = MAX( RR1, RR2 )
             ELSE if(MIN( RR1, RR2 )>0.0)then
                XNO3 = MIN( RR1, RR2 )
             ELSE!two negative roots !DS PW added 4th July 2012

                           !--------------------- return copied from above
                XNO3 = 0.0
                AH2O_High = 1000.0 * WH2O
                YNH4 = TWOSO4
                GNO3_High = HNO3
                ASO4_High = TSO4_HighA * MWSO4
                ANO3_High = NO3
                !ds ANH4 = YNH4 * MWNH4
                ANH4_High = min( YNH4 * MWNH4, TMASSNH3)  ! ds added "min"
                GNH3_High = TMASSNH3 - ANH4_High
                if( DEBUG%EQUIB .and. debug_flag ) WRITE( *, * ) ' TWO NEG ROOTS '
!                RETURN
                goto 333
               
             END IF
          ELSE
             XNO3 = - CC / BB   ! AA equals zero here
if( DEBUG%EQUIB .and. debug_flag ) print "(a,4es10.3)", "MARS NONDEGEN  ",  AA, BB, CC, XNO3
          END IF

          XNO3 = MIN( XNO3, TNO3 )

!...This version assumes no solid sulfate forms (supersaturated ) 
!...  Now update water

          CALL AWATER ( fRH, TSO4_HighA, YNH4, XNO3, AH2O)

!...ZSR relationship is used to set water levels. Units are
!...  10**(-6) kg water/ (cubic meter of air)
!...  The conversion from micromoles to moles is done by the units of WH2O.

          WH2O = 1.0E-3 * AH2O

!...Ionic balance determines the ammonium in solution.

          MAN = XNO3 / max(WH2O,FLOOR)
          MAS = TSO4_HighA / max(WH2O,FLOOR)
          MNH4 = 2.0 * MAS + MAN
          YNH4 = MNH4 * WH2O

 !st ...  FIXING
   if(MNH4<0. .or. MAS<0. .or. MAN<0.) then
      MNH4 = 1.e-30
      MAS  = 1.e-30
      MAN  = 1.e-30
   end if

!...MAS, MAN and MNH4 are the aqueous concentrations of sulfate, nitrate,
!...  and ammonium in molal units (moles/(kg water) ).

          STION = 3.0 * MAS + MAN
          CAT( 1 ) = 0.0
          CAT( 2 ) = MNH4 
          AN ( 1 ) = MAS
          AN ( 2 ) = MAN
          AN ( 3 ) = 0.0
          CALL ACTCOF ( CAT, AN, GAMS, MOLNU, PHIBAR , ERRMARK,1,debug_flag)
          GAMAAN = GAMS( 2, 2 )

!...Use GAMAAN for convergence control

          EROR = ABS( GAMOLD - GAMAAN ) / GAMOLD
          GAMOLD = GAMAAN

!...Check to see if we have a solution

          IF ( EROR <= TOLER1 ) THEN 
!!!            WRITE( 11, * ) RH, STION, GAMS( 1, 1 ),GAMS( 1, 2 ), GAMS( 1, 3 ),
!!!     &      GAMS( 2, 1 ), GAMS( 2, 2 ), GAMS( 2, 3 ), PHIBAR
! 2/25/99 IJA

            ASO4_High = TSO4_HighA * MWSO4
            ANO3_High = min(XNO3 * MWNO3,TMASSHNO3)
            ANH4_High = min( YNH4 * MWNH4, TMASSNH3 )
            GNO3_High = TMASSHNO3  - ANO3_High
            GNH3_High = TMASSNH3   - ANH4_High
            AH2O_High = 1000.0 * WH2O

!            RETURN
            goto 333
          END IF


1501    CONTINUE

!...If after NITR iterations no solution is found, then:
! FSB retain the initial values of nitrate particle and vapor
! 2/25/99 IJA
        ASO4_High = TSO4_HighA * MWSO4
        ANO3_High = NO3
        XNO3 = NO3 / MWNO3
        YNH4 = TWOSO4
!        ANH4 = YNH4 * MWNH4
        ANH4_High = min( YNH4 * MWNH4, TMASSNH3 )  ! ds pw added "min"
        CALL AWATER ( fRH, TSO4_HighA, YNH4, XNO3, AH2O_High)
        GNO3_High = HNO3
        GNH3_High = TMASSNH3 - ANH4_High
!        RETURN
        goto 333

333 CONTINUE !finished DO_RATIO_High_2

!      ELSE
     ENDIF
      IF ( DO_RATIO_Low_2 ) THEN    ! NH4/SO4 <= 2
       
!......................................
!......... Low Ammonia Case ...........
!......................................
      
!...coded by Dr. Francis S. Binkowski 12/8/91.(4/26/95)
      ! modified 8/28/98

!...All cases covered by this logic
 
        WH2O = 0.0
        CALL AWATER ( fRH, TSO4_LowA, TNH4, TNO3, AH2O )
        WH2O = 1.0E-3 * AH2O
        ZH2O = AH2O

!...convert 10**(-6) kg water/(cubic meter of air) to micrograms of water
!...  per cubic meter of air (1000 g = 1 kg)
! 2/25/99 IJA 
        ASO4_Low = TSO4_LowA * MWSO4
        ANH4_Low = TNH4 * MWNH4
        !dsSAFE ANO3 = NO3
        ANO3_Low = min( NO3, TMASSHNO3 )
        GNO3_Low = TMASSHNO3 - ANO3_Low
        GNH3_Low = FLOOR
        AH2O_Low = 1.0E3 *WH2O

!...Check for zero water.      

        !ds IF ( WH2O == 0.0 ) RETURN
        IF ( abs(WH2O) < FLOOR ) goto 111!RETURN
        ZSO4 = TSO4_LowA / max(WH2O,FLOOR)

!...ZSO4 is the molality of total sulfate i.e. MSO4 + MHSO4      

!!!         IF ( ZSO4 > 11.0 ) THEN

!...do not solve for aerosol nitrate for total sulfate molality
!...  greater than 11.0 because the model parameters break down
!...  greater than  9.0 because the model parameters break down

        IF ( ZSO4 > 9.0 ) THEN   ! 18 June 97
           goto 111
          !RETURN
        END IF

!...First solve with activity coeffs of 1.0, then iterate.

        PHIOLD = 1.0
        GAMANA = 1.0
        GAMAS1 = 1.0
        GAMAS2 = 1.0
        GAMAAB = 1.0
        GAMOLD = 1.0

!...All ammonia is considered to be aerosol ammonium. 

        MNH4 = TNH4 / max(WH2O,FLOOR)

!...MNH4 is the molality of ammonium ion.

        YNH4 = TNH4
      
!...loop for iteration
 
        DO 1601 NNN = 1, 150
          if( DEBUG%EQUIB ) then
            if(NNN > MAXNNN2 ) MAXNNN2 = NNN
            !if( MAXNNN2 > 140) print *, "NNN2 ", NNN, TNO3, TNH4, TWOSO4
          end if
          NITR = NNN

!...set up equilibrium constants including activities
!...  solve the system for hplus first then sulfate & nitrate

!GEOS added:
            IF ( .NOT. (& 
               IS_SAFE_DIV( GAMAS2, GAMAS1*GAMAS1*GAMAS1, R4=.FALSE. )& 
                .AND. IS_SAFE_DIV( KNA, GAMANA*GAMANA, R4=.FALSE. ) & 
                ) ) THEN
               
               if( nMarsErrors < 500 ) then
                  nMarsErrors = nMarsErrors + 1
                  WRITE(6,"(a,2i7,f7.2,4es11.2,14es12.3)") &
                    'RPMARES: not safe to divide...exit?', &
                     NNN, nMarsErrors, fRH, TSO4_LowA, TNH4, TNO3, AH2O, &
                     GAMAS2, GAMAS1*GAMAS1*GAMAS1, KNA, GAMANA*GAMANA,SO4, HNO3, NO3, NH3, NH4, RH, TEMP
               end if
!DSMARS                call CheckStop("UNDER/OVERFLOW in rpmares (MARS)")
!could continue: remove the stop and the warning if necessary

               GOTO 1601
            ENDIF

          RK2SA = K2SA * GAMAS2 * GAMAS2 / ( GAMAS1 * GAMAS1 * GAMAS1 )
          RKNA = KNA / ( GAMANA * GAMANA )
          RKNWET = RKNA * WH2O       
          T21  = ZSO4 - MNH4
          T221 = ZSO4 + T21

!...set up coefficients for cubic       

          A2 = RK2SA + RKNWET - T21
          A1 = RK2SA * RKNWET - T21 * ( RK2SA + RKNWET )    &
             - RK2SA * ZSO4 - RKNA * TNO3
          A0 = - (T21 * RK2SA * RKNWET                      &
             + RK2SA * RKNWET * ZSO4 + RK2SA * RKNA * TNO3 )   
         
          CALL CUBIC ( A2, A1, A0, NR, CRUTES )
       
!...Code assumes the smallest positive root is in CRUTES(1)
 
          HPLUS = CRUTES( 1 )
!GEOS added:
          IF (HPLUS <= 0d0)   goto 1601         

          BAL = HPLUS **3 + A2 * HPLUS**2 + A1 * HPLUS + A0
          MSO4 = RK2SA * ZSO4 / ( HPLUS + RK2SA )   ! molality of sulfate ion
          MHSO4 = ZSO4 - MSO4                       ! molality of bisulfate ion

!GEOS added:
          MHSO4 = MAX( 1.0d-10, MHSO4 ) !GEOS added

          MNA = RKNA * TNO3 / ( HPLUS + RKNWET )    ! molality of nitrate ion
          MNA = MAX( 0.0, MNA )
          MNA = MIN( MNA, TNO3 / max(WH2O,FLOOR) )
          XNO3 = MNA * WH2O
          !ds ANO3 = MNA * WH2O * MWNO3
          ANO3_Low = min( TMASSHNO3, MNA * WH2O * MWNO3)
! 2/25/99 IJA
          GNO3_Low = TMASSHNO3 - ANO3_Low

!GEOS added:
          ASO4_Low = MSO4 * WH2O * MWSO4 !pw added after [rjp, 12/12/01]

          if( DEBUG%EQUIB ) then
             if (GNO3_Low < 0.0 ) call CheckStop("NNN2 GNO3 NEG")
          end if
        
!...Calculate ionic strength      

          STION = 0.5 * ( HPLUS + MNA + MNH4 + MHSO4 + 4.0 * MSO4 )
          
!...Update water

          CALL AWATER ( fRH, TSO4_LowA, YNH4, XNO3, AH2O )

!...Convert 10**(-6) kg water/(cubic meter of air) to micrograms of water
!...  per cubic meter of air (1000 g = 1 kg)                       

          WH2O = 1.0E-3 * AH2O 
          CAT( 1 ) = HPLUS
          CAT( 2 ) = MNH4
          AN ( 1 ) = MSO4
          AN ( 2 ) = MNA
          AN ( 3 ) = MHSO4

          CALL ACTCOF ( CAT, AN, GAMS, MOLNU, PHIBAR, ERRMARK,2,debug_flag)

          GAMANA = GAMS( 1, 2 )
          GAMAS1 = GAMS( 1, 1 )
          GAMAS2 = GAMS( 1, 3 )
          GAMAAN = GAMS( 2, 2 )

          !Added by PW 11/12/2014 following comments
          !http://wiki.seas.harvard.edu/geos-chem/index.php?title=Aerosol_thermodynamical_equilibrium&redirect=no#RPMARES
          if ( ( abs( GAMANA ) < FLOOR ) .or. ( abs( GAMAS1 ) < FLOOR ) ) THEN
             goto 1601
          end if
          GAMAHAT = ( GAMAS2 * GAMAS2 / ( GAMAAB * GAMAAB ) )
          BHAT = KHAT * GAMAHAT 
!!!          EROR = ABS ( ( PHIOLD - PHIBAR ) / PHIOLD )
!!!          PHIOLD = PHIBAR
          EROR = ABS ( GAMOLD - GAMAHAT ) / GAMOLD 
          GAMOLD = GAMAHAT   

!...write out molalities and activity coefficient
!...  and return with good solution
          AH2O_Low = 1.0E3 *WH2O

          IF ( EROR <= TOLER2 ) THEN
!!!            WRITE(12,*) RH, STION,HPLUS,ZSO4,MSO4,MHSO4,MNH4,MNA
!!!            WRITE(11,*) RH, STION, GAMS(1,1),GAMS(1,2),GAMS(1,3),
!!!     &                  GAMS(2,1),GAMS(2,2),GAMS(2,3), PHIBAR

           goto 111
            !RETURN
          END IF

1601    CONTINUE     

!...after NITR iterations, failure to solve the system, no ANO3
! 2/25/99 IJA
        ANH4_Low = TNH4 * MWNH4
        GNH3_Low = FLOOR
        GNO3_Low = HNO3
        ANO3_Low = NO3
        ASO4_Low = TSO4 * MWSO4    ! PW after [rjp, 12/17/01]
        
        CALL AWATER ( fRH, TSO4_LowA, TNH4, TNO3, AH2O_Low )      
!        RETURN
           goto 111
111        continue
   
      END IF   ! ratio .gt. 2.0

      ASO4 = ASO4_High*High_Factor + ASO4_Low*(1.0-High_Factor)
      ANO3 = ANO3_High*High_Factor + ANO3_Low*(1.0-High_Factor)
      AH2O = AH2O_High*High_Factor + AH2O_Low*(1.0-High_Factor)
      ANH4 = ANH4_High*High_Factor + ANH4_Low*(1.0-High_Factor)
      GNH3 = GNH3_High*High_Factor + GNH3_Low*(1.0-High_Factor)
      GNO3 = GNO3_High*High_Factor + GNO3_Low*(1.0-High_Factor)

      end subroutine rpmares ! end RPMares
!<------------------------------------------------------------------------------->

      subroutine cubic(a2,a1,a0,nr,crutes)

  !.. subroutine  to find the roots of a cubic equation / 3rd order polynomial
  !.. formulae can be found in numer. recip.  on page 145
  !..  kiran  developed  this version on 25/4/1990
  !..  Dr. Francis Binkowski modified the routine on 6/24/91, 8/7/97
  !--------------------------------------------------------------

      implicit none

      real, intent(in)     :: a2,a1,a0
      integer, intent(out) :: nr      
      real, intent(out)    :: crutes(3)
!.. local
      real ::  qq,rr,a2sq,theta, sqrt3, one3rd
      real ::  dum1,dum2,part1,part2,part3,rrsq,phi,yy1,yy2,yy3
      real ::  costh, sinth

!emep:  7 digits not enough!      data sqrt3/1.732050808/, one3rd/0.333333333/

      sqrt3=sqrt(3.0)
      one3rd=1.0/3.0

!=======

      a2sq=a2*a2
      qq=(a2sq-3.*a1)/9.
      rr=( a2*(2.*a2sq - 9.*a1) + 27.*a0 )/54.
! CASE 1 THREE real ROOTS or  CASE 2 ONLY ONE real ROOT
      dum1=qq*qq*qq 
      rrsq=rr*rr
      dum2=dum1 - rrsq
      
      if(dum2 >= 0.) then
! NOW WE HAVE THREE real ROOTS
         phi=sqrt(dum1)
         if(abs(phi) <= 1.e-20) then 
!           write(10,*) ' cubic phi small, phi = ',phi
            crutes(1) = 0.0 
            crutes(2) = 0.0
            crutes(3) = 0.0
            nr = 0            
           stop 
         end if
         theta=acos(rr/phi)/3.0
         costh = cos(theta)
         sinth = sin(theta)
! *** use trig identities to simplify the expressions 
! *** binkowski's modification
         part1=sqrt(qq)
         yy1=part1*costh
         yy2=yy1-a2/3.0
         yy3=sqrt3*part1*sinth
         crutes(3) = -2.0*yy1 - a2/3.0
         crutes(2) = yy2 + yy3
         crutes(1) = yy2 - yy3
! *** SET NEGATIVE ROOTS TO A LARGE POSITIVE VALUE
         if(crutes(1) <= 0.0) crutes(1) = 1.0e9
         if(crutes(2) <= 0.0) crutes(2) =1.0e9
         if(crutes(3) <= 0.0) crutes(3) = 1.0e9
! *** put smallest positive root in crutes(1)
         crutes(1)=min( crutes(1),crutes(2),crutes(3))
         nr=3
      else  ! dum IS NEGATIVE
!     NOW HERE WE HAVE ONLY ONE real ROOT
         part1=sqrt(rrsq-dum1)
         part2=abs(rr)
         part3=(part1+part2)**one3rd
         crutes(1) = -sign(1.0,rr) * ( part3 + (qq/part3) ) - a2/3. 
         crutes(2)=0.
         crutes(3)=0.
!IAREV02...ADDITIONAL CHECK on NEGATIVE ROOTS
! *** SET NEGATIVE ROOTS TO A LARGE POSITIVE VALUE
!IA ACTIONIA
      if(crutes(1) <= 0.0) THEN
         crutes(1) = 1.0e9
  !!    if(debug_flag ) write(6,*) 'WARNING: NEGATIVE ROOTS IN CUBIC', crutes(1)
  !!st       stop
      end if
      nr=1
      end if

   end subroutine cubic
    
!>-------------------------------------------------------------------------------<
!<------------------------------------------------------------------------------->

      subroutine actcof ( CAT, AN, GAMA, MOLNU, PHIMULT , ERRMARK, IA2, debug_flag)

!C-----------------------------------------------------------------------
!C
!C DESCRIPTION:
!C
!C  This subroutine computes the activity coefficients of (2NH4+,SO4--),
!C  (NH4+,NO3-),(2H+,SO4--),(H+,NO3-),AND (H+,HSO4-) in aqueous
!C  multicomponent solution, using Bromley's model and Pitzer's method.
!C
!C REFERENCES:
!C
!C   Bromley, L.A. (1973) Thermodynamic properties of strong electrolytes
!C     in aqueous solutions.  AIChE J. 19, 313-320.
!C
!C   Chan, C.K. R.C. Flagen, & J.H.  Seinfeld (1992) Water Activities of 
!C     NH4NO3 / (NH4)2SO4 solutions, Atmos. Environ. (26A): 1661-1673.
!C
!C   Clegg, S.L. & P. Brimblecombe (1988) Equilibrium partial pressures 
!C     of strong acids over saline solutions - I HNO3, 
!C     Atmos. Environ. (22): 91-100
!C
!C   Clegg, S.L. & P. Brimblecombe (1990) Equilibrium partial pressures
!C     and mean activity and osmotic coefficients of 0-100% nitric acid 
!C     as a function of temperature,   J. Phys. Chem (94): 5369 - 5380
!C
!C   Pilinis, C. and J.H. Seinfeld (1987) Continued development of a
!C     general equilibrium model for inorganic multicomponent atmospheric
!C     aerosols.  Atmos. Environ. 21(11), 2453-2466.
!C


!
!CC ARGUMENT DESCRIPTION:
!
!C     CAT(1) : conc. of H+    (moles/kg)
!C     CAT(2) : conc. of NH4+  (moles/kg)
!C     AN(1)  : conc. of SO4-- (moles/kg)
!C     AN(2)  : conc. of NO3-  (moles/kg)
!C     AN(3)  : conc. of HSO4- (moles/kg)
!C     GAMA(2,1)    : mean molal ionic activity coeff for (2NH4+,SO4--)
!C     GAMA(2,2)    :  "    "     "       "       "    "  (NH4+,NO3-)
!C     GAMA(2,3)    :  "    "     "       "       "    "  (NH4+. HSO4-)
!C     GAMA(1,1)    :  "    "     "       "       "    "  (2H+,SO4--)
!C     GAMA(1,2)    :  "    "     "       "       "    "  (H+,NO3-)
!C     GAMA(1,3)    :  "    "     "       "       "    "  (H+,HSO4-)
!C     MOLNU   : the total number of moles of all ions.
!C     PHIMULT : the multicomponent paractical osmotic coefficient.
!C
!C REVISION HISTORY:
!C      Who       When        Detailed description of changes
!C   ---------   --------  -------------------------------------------
!C   S.Roselle   7/26/89   Copied parts of routine BROMLY, and began this
!C                         new routine using a method described by Pilinis
!C                         and Seinfeld 1987, Atmos. Envirn. 21 pp2453-2466.
!C   S.Roselle   7/30/97   Modified for use in Models-3
!C   F.Binkowski 8/7/97    Modified coefficients BETA0, BETA1, CGAMA
!C
!-----------------------------------------------------------------------




!...........INCLUDES and their descriptions

!      INCLUDE SUBST_XSTAT     ! M3EXIT status codes

!...........ARGUMENTS and their descriptions

            
     real, intent(in)  :: cat(2)    &  ! cation conc in moles/kg (input)
                         ,an (3)    &   ! anion conc in moles/kg (input)
                         ,molnu     &  ! tot # moles of all ions
                         ,phimult      ! multicomponent paractical osmotic coef
     real, intent(out) :: gama(2,3)   ! mean molal ionic activity coefs
     logical, intent(in) :: debug_flag

!....................................................................

      INTEGER    XSTAT0       ! Normal, successful completion
      PARAMETER (XSTAT0 = 0)
      INTEGER    XSTAT1       ! File I/O error
      PARAMETER (XSTAT1 = 1)
      INTEGER    XSTAT2       ! Execution error
      PARAMETER (XSTAT2 = 2)
      INTEGER    XSTAT3       ! Special  error
      PARAMETER (XSTAT3 = 3)
      INTEGER ERRMARK
      INTEGER IA2
      CHARACTER(len=120) :: XMSG

!...........PARAMETERS and their descriptions:

      INTEGER      NCAT                 ! number of cations
      PARAMETER  ( NCAT = 2 )

      INTEGER      NAN                  ! number of anions
      PARAMETER  ( NAN = 3 )


!...........SCRATCH LOCAL VARIABLES and their descriptions:

      CHARACTER(len=16), save :: PNAME            ! driver program name

      INTEGER      IAN                  ! anion indX
      INTEGER      ICAT                 ! cation indX

      REAL         FGAMA                ! 
      REAL         I                    ! ionic strength 
      REAL         R                    ! 
      REAL         S                    ! 
      REAL         TA                   ! 
      REAL         TB                   ! 
      REAL         TC                   ! 
      REAL         TEXPV                ! 
      REAL         TRM                  ! 
      REAL         TWOI                 ! 2*ionic strength
      REAL         TWOSRI               ! 2*sqrt of ionic strength
      REAL         ZBAR                 ! 
      REAL         ZBAR2                ! 
      REAL         ZOT1                 ! 
      REAL         SRI                  ! square root of ionic strength 
      REAL         F2( NCAT )           ! 
      REAL         F1( NAN )            ! 
      REAL         ZP( NCAT )           ! absolute value of charges of cation
      REAL         ZM( NAN )            ! absolute value of charges of anion
      REAL         BGAMA ( NCAT, NAN )  ! 
      REAL         X     ( NCAT, NAN )  ! 
      REAL         M     ( NCAT, NAN )  ! molality of each electrolyte
      REAL         LGAMA0( NCAT, NAN )  ! binary activity coefficients
      REAL         Y     ( NAN, NCAT )  ! 
      REAL         BETA0 ( NCAT, NAN )  ! binary activity coefficient parameter
      REAL         BETA1 ( NCAT, NAN )  ! binary activity coefficient parameter
      REAL         CGAMA ( NCAT, NAN )  ! binary activity coefficient parameter
      REAL         V1    ( NCAT, NAN )  ! number of cations in electrolyte formula
      REAL         V2    ( NCAT, NAN )  ! number of anions in electrolyte formula

      DATA         ZP / 1.0, 1.0 /
      DATA         ZM / 2.0, 1.0, 1.0 /
      DATA         XMSG / ' ' /
      DATA         PNAME / 'ACTCOF' /

! *** Sources for the coefficients BETA0, BETA1, CGAMA:
 
! *** (1,1);(1,3)  - Clegg & Brimblecombe (1988)
! *** (2,3)        - Pilinis & Seinfeld (1987), cgama different 
! *** (1,2)        - Clegg & Brimblecombe (1990)
! *** (2,1);(2,2)  - Chan, Flagen & Seinfeld (1992)                                 
      
! *** now set the basic constants, BETA0, BETA1, CGAMA 

      DATA BETA0(1,1) /2.98E-2/,    BETA1(1,1) / 0.0/,          &
           CGAMA(1,1) / 4.38E-2/                                ! 2H+SO4-
     
      DATA BETA0(1,2) / 1.2556E-1/,   BETA1(1,2) / 2.8778E-1/,  & 
           CGAMA(1,2) / -5.59E-3/                               ! HNO3
     
      DATA BETA0(1,3) / 2.0651E-1/,   BETA1(1,3) / 5.556E-1/,   &
           CGAMA(1,3) /0.0/                                     ! H+HSO4-
     
      DATA BETA0(2,1) /4.6465E-2/,   BETA1(2,1) /-0.54196/,     &    
           CGAMA(2,1) /-1.2683E-3/                              ! (NH4)2SO4
     
      DATA BETA0(2,2) /-7.26224E-3/, BETA1(2,2) /-1.168858/,    &
           CGAMA(2,2) /3.51217E-5/                              ! NH4NO3
     
      DATA BETA0(2,3) / 4.494E-2/,    BETA1(2,3) / 2.3594E-1/,  &
           CGAMA(2,3) /-2.962E-3/                               ! NH4HSO4

      DATA V1(1,1), V2(1,1) / 2.0, 1.0 /     ! 2H+SO4-
      DATA V1(2,1), V2(2,1) / 2.0, 1.0 /     ! (NH4)2SO4
      DATA V1(1,2), V2(1,2) / 1.0, 1.0 /     ! HNO3 
      DATA V1(2,2), V2(2,2) / 1.0, 1.0 /     ! NH4NO3
      DATA V1(1,3), V2(1,3) / 1.0, 1.0 /     ! H+HSO4-
      DATA V1(2,3), V2(2,3) / 1.0, 1.0 /     ! NH4HSO4

!-----------------------------------------------------------------------
!  begin body of subroutine ACTCOF

!...compute ionic strength
      I = 0.0

      DO ICAT = 1, NCAT
        I = I + CAT( ICAT ) * ZP( ICAT ) * ZP( ICAT )
      END DO

      DO IAN = 1, NAN
        I = I + AN( IAN ) * ZM( IAN ) * ZM( IAN )
      END DO

      I = 0.5 * I

!...check for problems in the ionic strength

      IF ( I .EQ. 0.0 ) THEN

        DO IAN = 1, NAN
          DO ICAT = 1, NCAT
            GAMA( ICAT, IAN ) = 0.0
          END DO
        END DO
        
        XMSG = 'Ionic strength is zero...returning zero activities'
       if(debug_flag ) WRITE(6,*) XMSG 
        RETURN

      ELSE IF ( I .LT. 0.0 ) THEN
        XMSG = 'Ionic strength below zero...negative concentrations'
   if(debug_flag ) then
        WRITE(6,*) XMSG
        WRITE(6,*) 'called over ', IA2
        WRITE(6,*) ' I =', I
        WRITE(6,*) 'CAT=', CAT
        WRITE(6,*) 'AN=', AN
        WRITE(6,*) 'GAMA=', GAMA
        WRITE(6,*) 'MOLNU=',MOLNU
        WRITE(6,*) 'PHIMULT=',PHIMULT
   end if
 !!       CALL M3EXIT( PNAME, 0, 0, XMSG, XSTAT2 )
 !emep1.2      call stop_test(.true.,me,NPROC,ios,'##MARS-negat.con')
     END IF

!...compute some essential expressions

      SRI    = SQRT( I )
      TWOSRI = 2.0 * SRI
      TWOI   = 2.0 * I
      TEXPV  = 1.0 - EXP( -TWOSRI ) * ( 1.0 + TWOSRI - TWOI )
      R      = 1.0 + 0.75 * I
      S      = 1.0 + 1.5  * I
      ZOT1   = 0.511 * SRI / ( 1.0 + SRI )

!...Compute binary activity coeffs

      FGAMA = -0.392 * ( ( SRI / ( 1.0 + 1.2 * SRI )         &
            + ( 2.0 / 1.2 ) * ALOG( 1.0 + 1.2 * SRI ) ) )

      DO ICAT = 1, NCAT
        DO IAN = 1, NAN

         BGAMA( ICAT, IAN ) = 2.0 * BETA0( ICAT, IAN )        &       
                            + ( 2.0 * BETA1( ICAT, IAN ) / ( 4.0 * I ) )   &
                            * TEXPV

!...compute the molality of each electrolyte for given ionic strength

          M( ICAT, IAN ) = ( CAT( ICAT )**V1( ICAT, IAN )       &
                         *   AN( IAN )**V2( ICAT, IAN ) )**( 1.0     &
                         / ( V1( ICAT, IAN ) + V2( ICAT, IAN ) ) )   

!...calculate the binary activity coefficients

         LGAMA0( ICAT, IAN ) = ( ZP( ICAT ) * ZM( IAN ) * FGAMA       &
                             + M( ICAT, IAN )                         &
                             * ( 2.0 * V1( ICAT, IAN ) * V2( ICAT, IAN )   &
                             / ( V1( ICAT, IAN ) + V2( ICAT, IAN ) )       &
                             * BGAMA( ICAT, IAN ) )                    &
                             + M( ICAT, IAN ) * M( ICAT, IAN )         &
                             * ( 2.0 * ( V1( ICAT, IAN )               &
                             * V2( ICAT, IAN ) )**1.5                  &
                             / ( V1( ICAT, IAN ) + V2( ICAT, IAN ) )    &
                             * CGAMA( ICAT, IAN ) ) ) / 2.302585093

        END DO
      END DO

!...prepare variables for computing the multicomponent activity coeffs

      DO IAN = 1, NAN
        DO ICAT = 1, NCAT
          ZBAR = ( ZP( ICAT ) + ZM( IAN ) ) * 0.5
          ZBAR2 = ZBAR * ZBAR
          Y( IAN, ICAT ) = ZBAR2 * AN( IAN ) / I
          X( ICAT, IAN ) = ZBAR2 * CAT( ICAT ) / I
        END DO
      END DO

      DO IAN = 1, NAN
        F1( IAN ) = 0.0
        DO ICAT = 1, NCAT
          F1( IAN ) = F1( IAN ) + X( ICAT, IAN ) * LGAMA0( ICAT, IAN )    &
                    + ZOT1 * ZP( ICAT ) * ZM( IAN ) * X( ICAT, IAN )
        END DO
      END DO

      DO ICAT = 1, NCAT
        F2( ICAT ) = 0.0
        DO IAN = 1, NAN
          F2( ICAT ) = F2( ICAT ) + Y( IAN, ICAT ) * LGAMA0( ICAT, IAN )    &
                     + ZOT1 * ZP( ICAT ) * ZM( IAN ) * Y( IAN, ICAT )
        END DO
      END DO

!...now calculate the multicomponent activity coefficients

      DO IAN = 1, NAN
        DO ICAT = 1, NCAT

          TA = -ZOT1 * ZP( ICAT ) * ZM( IAN )
          TB = ZP( ICAT ) * ZM( IAN ) / ( ZP( ICAT ) + ZM( IAN ) )
          TC = ( F2( ICAT ) / ZP( ICAT ) + F1( IAN ) / ZM( IAN ) )
          TRM = TA + TB * TC

          IF ( TRM > 30.0 ) THEN
            GAMA( ICAT, IAN ) = 1.0E+30
            XMSG = 'Multicomponent activity coefficient is >>'
       !!     if(debug_flag )  WRITE(6,*) XMSG, gama(icat,ian)
            ERRMARK=2
 
          ELSE
            GAMA( ICAT, IAN ) = 10.0**TRM
          END IF

        END DO
      END DO
 
    end subroutine actcof 
!>-------------------------------------------------------------------------------<

!routine to interface RPMARES_new with emep model
    SUBROUTINE DO_RPMARES_new( SO4, HNO3, NO3, NH3, NH4, RH, TEMP,   &
         ASO4, ANO3, AH2O, ANH4, GNH3, GNO3,   &
         ERRMARK,debug_flag) 
      INTEGER     ERRMARK
      logical, intent(in) :: debug_flag
      
  real             ::  SO4   &     ! Total sulfate in micrograms / m**3 
                      ,HNO3  &     ! Total nitric acid in micrograms / m**3
                      ,NO3   &     ! Total nitrate in micrograms / m**3
                      ,NH3   &     ! Total ammonia in micrograms / m**3
                      ,NH4   &     ! Total ammonium in micrograms / m**3
                      ,RH    &     ! Fractional relative humidity 
                      ,TEMP        ! Temperature in Kelvin 

  real             ::  ASO4  &     ! Aerosol sulfate in micrograms / m**3 
                      ,ANO3  &     ! Aerosol nitrate in micrograms / m**3
                      ,AH2O  &     ! Aerosol liquid water content water in micrograms / m**3
                      ,ANH4  &     ! Aerosol ammonium in micrograms / m**3
                      ,GNO3  &     ! Gas-phase nitric acid in micrograms / m**3
                      ,GNH3        ! Gas-phase ammonia in micrograms / m**3

!not used
real  AHSO4 ! Aerosol phase in bisulfate in MICROGRAMS/M**3 
!

       ANH4 = NH4
       GNH3 = NH3
       ANO3 = NO3
       GNO3 = HNO3
      
       call RPMARES_new( SO4,  GNO3,  GNH3, RH,   TEMP,&
                         ASO4, AHSO4, ANO3, AH2O, ANH4 )


     end SUBROUTINE DO_RPMARES_new
!from GEOS-Chem

      SUBROUTINE RPMARES_new( SO4,  GNO3,  GNH3, RH,   TEMP,&
                         ASO4, AHSO4, ANO3, AH2O, ANH4 )
!
!******************************************************************************
!
! Description:
!
!   ARES calculates the chemical composition of a sulfate/nitrate/
!   ammonium/water aerosol based on equilibrium thermodynamics.
!
!   This code considers two regimes depending upon the molar ratio
!   of ammonium to sulfate.
!
!   For values of this ratio less than 2,the code solves a cubic for
!   hydrogen ion molality, H+,  and if enough ammonium and liquid
!   water are present calculates the dissolved nitric acid. For molal
!   ionic strengths greater than 50, nitrate is assumed not to be present.
!
!   For values of the molar ratio of 2 or greater, all sulfate is assumed
!   to be ammonium sulfate and a calculation is made for the presence of
!   ammonium nitrate.
!
!   The Pitzer multicomponent approach is used in subroutine ACTCOF to
!   obtain the activity coefficients. Abandoned -7/30/97 FSB
!
!   The Bromley method of calculating the multicomponent activity coefficients
!    is used in this version 7/30/97 SJR/FSB
!
!   The calculation of liquid water
!   is done in subroutine water. Details for both calculations are given
!   in the respective subroutines.
!
!   Based upon MARS due to
!   P. Saxena, A.B. Hudischewskyj, C. Seigneur, and J.H. Seinfeld,
!   Atmos. Environ., vol. 20, Number 7, Pages 1471-1483, 1986.
!
!   and SCAPE due to
!   Kim, Seinfeld, and Saxeena, Aerosol Sience and Technology,
!   Vol 19, number 2, pages 157-181 and pages 182-198, 1993.
!
! NOTE: All concentrations supplied to this subroutine are TOTAL
!       over gas and aerosol phases
!
! Parameters:
!
!  SO4   : Total sulfate in MICROGRAMS/M**3 as sulfate 
!  GNO3  : Nitric Acid vapor in MICROGRAMS/M**3 as nitric acid 
!  GNH3  : Gas phase ammonia in MICROGRAMS/M**3 
!  RH    : Fractional relative humidity 
!  TEMP  : Temperature in Kelvin 
!  ASO4  : Aerosol phase sulfate in MICROGRAMS/M**3 
!  AHSO4 : Aerosol phase in bisulfate in MICROGRAMS/M**3 [rjp, 12/12/01]
!  ANO3  : Aerosol phase nitrate in MICROGRAMS/M**3 
!  ANH4  : Aerosol phase ammonium in MICROGRAMS/M**3 
!  AH2O  : Aerosol phase water in MICROGRAMS/M**3 
!
! Revision History:
!      Who       When        Detailed description of changes
!   ---------   --------  -------------------------------------------
!   S.Roselle   11/10/87  Received the first version of the MARS code
!   S.Roselle   12/30/87  Restructured code
!   S.Roselle   2/12/88   Made correction to compute liquid-phase
!                         concentration of H2O2.
!   S.Roselle   5/26/88   Made correction as advised by SAI, for
!                         computing H+ concentration.
!   S.Roselle   3/1/89    Modified to operate with EM2
!   S.Roselle   5/19/89   Changed the maximum ionic strength from
!                         100 to 20, for numerical stability.
!   F.Binkowski 3/3/91    Incorporate new method for ammonia rich case
!                         using equations for nitrate budget.
!   F.Binkowski 6/18/91   New ammonia poor case which
!                         omits letovicite.
!   F.Binkowski 7/25/91   Rearranged entire code, restructured
!                         ammonia poor case.
!   F.Binkowski 9/9/91    Reconciled all cases of ASO4 to be output
!                         as SO4--
!   F.Binkowski 12/6/91   Changed the ammonia defficient case so that
!                         there is only neutralized sulfate (ammonium
!                         sulfate) and sulfuric acid.
!   F.Binkowski 3/5/92    Set RH bound on AWAS to 37 % to be in agreement
!                          with the Cohen et al. (1987)  maximum molality
!                          of 36.2 in Table III.( J. Phys Chem (91) page
!                          4569, and Table IV p 4587.)
!   F.Binkowski 3/9/92    Redid logic for ammonia defficient case to remove
!                         possibility for denomenator becoming zero;
!                         this involved solving for H+ first.
!                         Note that for a relative humidity
!                          less than 50%, the model assumes that there is no
!                          aerosol nitrate.
!   F.Binkowski 4/17/95   Code renamed  ARES (AeRosol Equilibrium System)
!                          Redid logic as follows
!                         1. Water algorithm now follows Spann & Richardson
!                         2. Pitzer Multicomponent method used
!                         3. Multicomponent practical osmotic coefficient
!                            use to close iterations.
!                         4. The model now assumes that for a water
!                            mass fraction WFRAC less than 50% there is
!                            no aerosol nitrate.
!   F.Binkowski 7/20/95   Changed how nitrate is calculated in ammonia poor
!                         case, and changed the WFRAC criterion to 40%.
!                         For ammonium to sulfate ratio less than 1.0
!                         all ammonium is aerosol and no nitrate aerosol
!                         exists.
!   F.Binkowski 7/21/95   Changed ammonia-ammonium in ammonia poor case to
!                         allow gas-phase ammonia to exist.
!   F.Binkowski 7/26/95   Changed equilibrium constants to values from
!                         Kim et al. (1993)
!   F.Binkowski 6/27/96   Changed to new water format
!   F.Binkowski 7/30/97   Changed to Bromley method for multicomponent
!                         activity coefficients. The binary activity
!                         coefficients
!                         are the same as the previous version
!   F.Binkowski 8/1/97    Changed minimum sulfate from 0.0 to 1.0e-6 i.e.
!                         1 picogram per cubic meter
!   F.Binkowski 2/23/98   Changes to code made by Ingmar Ackermann to
!                         deal with precision problems on workstations 
!                         incorporated in to this version.  Also included
!                         are his improved descriptions of variables. 
!  F. Binkowski 8/28/98   changed logic as follows: 
!                         If iterations fail, initial values of nitrate
!                          are retained. 
!                         Total mass budgets are changed to account for gas
!                         phase returned.
!  F.Binkowski 10/01/98   Removed setting RATIO to 5 for low to 
!                         to zero sulfate sulfate case.
!  F.Binkowski 01/10/2000 reconcile versions
!
!  F.Binkowski 05/17/2000 change to logic for calculating RATIO
!  F.Binkowski 04/09/2001 change for very low values of RATIO,
!                         RATIO < 0.5, no iterative calculations are done
!                         in low ammonia case a MAX(1.0e-10, MSO4) IS
!                         applied, and the iteration count is
!                         reduced to fifty for each iteration loop.
!  R. Yantosca 09/25/2002 Bundled into "rpmares_mod.f".  Declared all REALs
!                          as REAL*8's.  Cleaned up comments.  Also now force
!                          double precision explicitly with "D" exponents.
!  P. Le Sager and        Bug fix for low ammonia case -- prevent floating
!  R. Yantosca 04/10/2008  point underflow and NaN's.
!  P. Le Sager 06/10/2008 Better catch of over/underflow for low ammonia case
!******************************************************************************
!
      ! References to F90 modules
!      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP, IS_SAFE_DIV

      !=================================================================
      ! ARGUMENTS and their descriptions
      !=================================================================
      REAL*8 :: SO4              ! Total sulfate in micrograms / m**3
      REAL*8 :: GNO3             ! Gas-phase nitric acid in micrograms / m**3
      REAL*8 :: GNH3             ! Gas-phase ammonia in micrograms / m**3 
      REAL*8 :: RH               ! Fractional relative humidity
      REAL*8 :: TEMP             ! Temperature in Kelvin
      REAL*8 :: ASO4             ! Aerosol sulfate in micrograms / m**3
      REAL*8 :: AHSO4            ! Aerosol bisulfate in micrograms / m**3
      REAL*8 :: ANO3             ! Aerosol nitrate in micrograms / m**3
      REAL*8 :: AH2O             ! Aerosol liquid water content water in
                                 !   micrograms / m**3
      REAL*8 :: ANH4             ! Aerosol ammonium in micrograms / m**3

      !=================================================================
      ! PARAMETERS and their descriptions:
      !=================================================================

      ! Molecular weights
      REAL*8, PARAMETER :: MWNACL = 58.44277d0               ! NaCl
      REAL*8, PARAMETER :: MWNO3  = 62.0049d0                ! NO3
      REAL*8, PARAMETER :: MWHNO3 = 63.01287d0               ! HNO3
      REAL*8, PARAMETER :: MWSO4  = 96.0576d0                ! SO4
      REAL*8, PARAMETER :: MWHSO4 = MWSO4 + 1.0080d0         ! HSO4
      REAL*8, PARAMETER :: MH2SO4 = 98.07354d0               ! H2SO4
      REAL*8, PARAMETER :: MWNH3  = 17.03061d0               ! NH3
      REAL*8, PARAMETER :: MWNH4  = 18.03858d0               ! NH4
      REAL*8, PARAMETER :: MWORG  = 16.0d0                   ! Organic Species
      REAL*8, PARAMETER :: MWCL   = 35.453d0                 ! Chloride
      REAL*8, PARAMETER :: MWAIR  = 28.964d0                 ! AIR
      REAL*8, PARAMETER :: MWLCT  = 3.0d0 * MWNH4 +   &       ! Letovicite
                                   2.0d0 * MWSO4 + 1.0080d0  
      REAL*8, PARAMETER :: MWAS   = 2.0d0 * MWNH4 + MWSO4    ! Amm. Sulfate
      REAL*8, PARAMETER :: MWABS  = MWNH4 + MWSO4 + 1.0080d0 ! Amm. Bisulfate

      ! Minimum value of sulfate aerosol concentration
      REAL*8, PARAMETER :: MINSO4 = 1.0d-6 / MWSO4 

      ! Minimum total nitrate cncentration
      REAL*8, PARAMETER :: MINNO3 = 1.0d-6 / MWNO3  

      ! Force a minimum concentration
      REAL*8, PARAMETER :: FLOOR  = 1.0d-30 

      ! Tolerances for convergence test.  NOTE: We now have made these
      ! parameters so they don't lose their values (phs, bmy, 4/10/08)
      REAL*8, PARAMETER :: TOLER1 = 0.00001d0       
      REAL*8, PARAMETER :: TOLER2 = 0.001d0       

      !=================================================================
      ! SCRATCH LOCAL VARIABLES and their descriptions:
      !=================================================================

      INTEGER :: IRH              ! Index set to percent relative humidity
      INTEGER :: NITR             ! Number of iterations for activity
                                  !   coefficients
      INTEGER :: NNN              ! Loop index for iterations
      INTEGER :: NR               ! Number of roots to cubic equation for
                                  ! H+ ciaprecision
      REAL*8  :: A0               ! Coefficients and roots of
      REAL*8  :: A1               ! Coefficients and roots of
      REAL*8  :: A2               ! Coefficients and roots of
      REAL*8  :: AA               ! Coefficients and discriminant for
                                  ! quadratic equation for ammonium nitrate
      REAL*8  :: BAL              ! internal variables ( high ammonia case)
      REAL*8  :: BB               ! Coefficients and discriminant for
                                  !   quadratic equation for ammonium nitrate
      REAL*8  :: BHAT             ! Variables used for ammonia solubility
      REAL*8  :: CC               ! Coefficients and discriminant for
                                  !   quadratic equation for ammonium nitrate
      REAL*8  :: CONVT            ! Factor for conversion of units
      REAL*8  :: DD               ! Coefficients and discriminant for
                                  !   quadratic equation for ammonium nitrate
      REAL*8  :: DISC             ! Coefficients and discriminant for
                                  !   quadratic equation for ammonium nitrate
      REAL*8  :: EROR             ! Relative error used for convergence test
      REAL*8  :: FNH3             ! "Free ammonia concentration", that
                                  !   which exceeds TWOSO4
      REAL*8  :: GAMAAB           ! Activity Coefficient for (NH4+,
                                  !   HSO4-)GAMS( 2,3 )
      REAL*8  :: GAMAAN           ! Activity coefficient for (NH4+, NO3-)
                                  !   GAMS( 2,2 )
      REAL*8  :: GAMAHAT          ! Variables used for ammonia solubility
      REAL*8  :: GAMANA           ! Activity coefficient for (H+ ,NO3-)
                                  !   GAMS( 1,2 )
      REAL*8  :: GAMAS1           ! Activity coefficient for (2H+, SO4--)
                                  !   GAMS( 1,1 )
      REAL*8  :: GAMAS2           ! Activity coefficient for (H+, HSO4-)
                                  !   GAMS( 1,3 )
      REAL*8  :: GAMOLD           ! used for convergence of iteration
      REAL*8  :: GASQD            ! internal variables ( high ammonia case)
      REAL*8  :: HPLUS            ! Hydrogen ion (low ammonia case) (moles
                                  !   / kg water)
      REAL*8  :: K1A              ! Equilibrium constant for ammonia to
                                  !   ammonium
      REAL*8  :: K2SA             ! Equilibrium constant for
                                  !   sulfate-bisulfate (aqueous)
      REAL*8  :: K3               ! Dissociation constant for ammonium
                                  !   nitrate
      REAL*8  :: KAN              ! Equilibrium constant for ammonium
                                  !   nitrate (aqueous)
      REAL*8  :: KHAT             ! Variables used for ammonia solubility
      REAL*8  :: KNA              ! Equilibrium constant for nitric acid
                                  !   (aqueous)
      REAL*8  :: KPH              ! Henry's Law Constant for ammonia
      REAL*8  :: KW               ! Equilibrium constant for water
                                  !  dissociation
      REAL*8  :: KW2              ! Internal variable using KAN
      REAL*8  :: MAN              ! Nitrate (high ammonia case) (moles /
                                  !   kg water)
      REAL*8  :: MAS              ! Sulfate (high ammonia case) (moles /
                                  !   kg water)
      REAL*8  :: MHSO4            ! Bisulfate (low ammonia case) (moles /
                                  !   kg water)
      REAL*8  :: MNA              ! Nitrate (low ammonia case) (moles / kg
                                  !   water)
      REAL*8  :: MNH4             ! Ammonium (moles / kg water)
      REAL*8  :: MOLNU            ! Total number of moles of all ions
      REAL*8  :: MSO4             ! Sulfate (low ammonia case) (moles / kg
                                  !   water)
      REAL*8  :: PHIBAR           ! Practical osmotic coefficient
      REAL*8  :: PHIOLD           ! Previous value of practical osmotic
                                  !   coefficient used for convergence of
                                  !   iteration
      REAL*8  :: RATIO            ! Molar ratio of ammonium to sulfate
      REAL*8  :: RK2SA            ! Internal variable using K2SA
      REAL*8  :: RKNA             ! Internal variables using KNA
      REAL*8  :: RKNWET           ! Internal variables using KNA
      REAL*8  :: RR1
      REAL*8  :: RR2
      REAL*8  :: STION            ! Ionic strength
      REAL*8  :: T1               ! Internal variables for temperature
                                  !   corrections
      REAL*8  :: T2               ! Internal variables for temperature
                                  !   corrections
      REAL*8  :: T21              ! Internal variables of convenience (low
                                  !   ammonia case)
      REAL*8  :: T221             ! Internal variables of convenience (low
                                  !   ammonia case)
      REAL*8  :: T3               ! Internal variables for temperature
                                  !   corrections
      REAL*8  :: T4               ! Internal variables for temperature
                                  !   corrections
      REAL*8  :: T6               ! Internal variables for temperature
                                  !   corrections
      REAL*8  :: TNH4             ! Total ammonia and ammonium in
                                  !   micromoles / meter ** 3
      REAL*8  :: TNO3             ! Total nitrate in micromoles / meter ** 3
      REAL*8  :: TSO4             ! Total sulfate in micromoles / meter ** 3
      REAL*8  :: TWOSO4           ! 2.0 * TSO4  (high ammonia case) (moles
                                  !   / kg water)
      REAL*8  :: WFRAC            ! Water mass fraction
      REAL*8  :: WH2O             ! Aerosol liquid water content (internally)
                                  !   micrograms / meter **3 on output
                                  !   internally it is 10 ** (-6) kg (water)
                                  !   / meter ** 3
                                  !   the conversion factor (1000 g = 1 kg)
                                  !   is applied for AH2O output
      REAL*8  :: WSQD             ! internal variables ( high ammonia case)
      REAL*8  :: XNO3             ! Nitrate aerosol concentration in
                                  ! micromoles / meter ** 3
      REAL*8  :: XXQ              ! Variable used in quadratic solution
      REAL*8  :: YNH4             ! Ammonium aerosol concentration in
                                  !  micromoles / meter** 3
      REAL*8  :: ZH2O             ! Water variable saved in case ionic
                                  !  strength too high.
      REAL*8  :: ZSO4             ! Total sulfate molality - mso4 + mhso4
                                  !  (low ammonia case) (moles / kg water)
      REAL*8  :: CAT( 2 )         ! Array for cations (1, H+); (2, NH4+)
                                  !  (moles / kg water)
      REAL*8  :: AN ( 3 )         ! Array for anions (1, SO4--); (2,
                                  !   NO3-); (3, HSO4-)  (moles / kg water)
      REAL*8  :: CRUTES( 3 )      ! Coefficients and roots of
      REAL*8  :: GAMS( 2, 3 )     ! Array of activity coefficients
      REAL*8  :: TMASSHNO3        ! Total nitrate (vapor and particle) 
      REAL*8  :: GNO3_IN, ANO3_IN                
      
      !=================================================================
      ! RPMARES begins here!
      ! convert into micromoles/m**3
      !=================================================================

      ! For extremely low relative humidity ( less than 1% ) set the 
      ! water content to a minimum and skip the calculation.
      IF ( RH .LT. 0.01 ) THEN
         AH2O = FLOOR
         RETURN
      ENDIF 

      ! total sulfate concentration
      TSO4 = MAX( FLOOR, SO4 / MWSO4  )      
      ASO4 = SO4

      !Cia models3 merge NH3/NH4 , HNO3,NO3 here
      !c *** recommended by Dr. Ingmar Ackermann

      ! total nitrate
      TNO3      = MAX( 0.0d0, ( ANO3 / MWNO3 + GNO3 / MWHNO3 ) )            

      ! total ammonia
      TNH4      = MAX( 0.0d0, ( GNH3 / MWNH3 + ANH4 / MWNH4 )  )

      GNO3_IN   = GNO3
      ANO3_IN   = ANO3
      TMASSHNO3 = MAX( 0.0d0, GNO3 + ANO3 )     

      ! set the  molar ratio of ammonium to sulfate
      RATIO = TNH4 / TSO4

      ! validity check for negative concentration
      IF ( TSO4 < 0.0d0 .OR. TNO3 < 0.0d0 .OR. TNH4 < 0.0d0 ) THEN
          PRINT*, 'TSO4 : ', TSO4
          PRINT*, 'TNO3 : ', TNO3
          PRINT*, 'TNH4 : ', TNH4
          !CALL GEOS_CHEM_STOP
          stop
      ENDIF   

      ! now set humidity index IRH as a percent
      IRH = NINT( 100.0 * RH )

      ! now set humidity index IRH as a percent
      IRH = MAX(  1, IRH )
      IRH = MIN( 99, IRH )

      !=================================================================
      ! Specify the equilibrium constants at  correct temperature.  
      ! Also change units from ATM to MICROMOLE/M**3 (for KAN, KPH, and K3 )
      ! Values from Kim et al. (1993) except as noted.
      ! Equilibrium constant in Kim et al. (1993)
      !   K = K0 exp[ a(T0/T -1) + b(1+log(T0/T)-T0/T) ], T0 = 298.15 K
      !   K = K0 EXP[ a T3 + b T4 ] in the code here.
      !=================================================================
      CONVT = 1.0d0 / ( 0.082d0 * TEMP )
      T6    = 0.082d-9 *  TEMP
      T1    = 298.0d0 / TEMP
      T2    = LOG( T1 )
      T3    = T1 - 1.0d0
      T4    = 1.0d0 + T2 - T1

      !=================================================================
      ! Equilibrium Relation
      ! 
      ! HSO4-(aq)         = H+(aq)   + SO4--(aq)  ; K2SA
      ! NH3(g)            = NH3(aq)               ; KPH
      ! NH3(aq) + H2O(aq) = NH4+(aq) + OH-(aq)    ; K1A
      ! HNO3(g)           = H+(aq)   + NO3-(aq)   ; KNA
      ! NH3(g) + HNO3(g)  = NH4NO3(s)             ; K3
      ! H2O(aq)           = H+(aq)   + OH-(aq)    ; KW
      !=================================================================
      KNA  = 2.511d+06 *  EXP(  29.17d0 * T3 + 16.83d0 * T4 ) * T6
      K1A  = 1.805d-05 *  EXP(  -1.50d0 * T3 + 26.92d0 * T4 )
      K2SA = 1.015d-02 *  EXP(   8.85d0 * T3 + 25.14d0 * T4 )
      KW   = 1.010d-14 *  EXP( -22.52d0 * T3 + 26.92d0 * T4 )
      KPH  = 57.639d0  *  EXP(  13.79d0 * T3 -  5.39d0 * T4 ) * T6
      !K3   =  5.746E-17 * EXP( -74.38 * T3 + 6.12  * T4 ) * T6 * T6
      KHAT =  KPH * K1A / KW
      KAN  =  KNA * KHAT

      ! Compute temperature dependent equilibrium constant for NH4NO3
      ! (from Mozurkewich, 1993)
      K3 = EXP( 118.87d0  - 24084.0d0 / TEMP -  6.025d0  * LOG( TEMP ) )

      ! Convert to (micromoles/m**3) **2
      K3     = K3 * CONVT * CONVT

      WH2O   = 0.0d0
      STION  = 0.0d0
      AH2O   = 0.0d0
      MAS    = 0.0d0
      MAN    = 0.0d0
      HPLUS  = 0.0d0
      NITR   = 0
      NR     = 0
      GAMAAN = 1.0d0
      GAMOLD = 1.0d0

      ! If there is very little sulfate and  nitrate 
      ! set concentrations to a very small value and return.
      IF ( ( TSO4 .LT. MINSO4 ) .AND. ( TNO3 .LT. MINNO3 ) ) THEN
         ASO4  = MAX( FLOOR, ASO4  )
         AHSO4 = MAX( FLOOR, AHSO4 ) ! [rjp, 12/12/01]
         ANO3  = MAX( FLOOR, ANO3  )
         ANH4  = MAX( FLOOR, ANH4  )
         WH2O  = FLOOR
         AH2O  = FLOOR
         GNH3  = MAX( FLOOR, GNH3  )
         GNO3  = MAX( FLOOR, GNO3  )
         
         RETURN
      ENDIF

      !=================================================================
      ! High Ammonia Case
      !=================================================================
      IF ( RATIO .GT. 2.0d0 ) THEN
        
         GAMAAN = 0.1d0

         ! Set up twice the sulfate for future use.
         TWOSO4 = 2.0d0 * TSO4
         XNO3   = 0.0d0
         YNH4   = TWOSO4

         ! Treat different regimes of relative humidity
         !
         ! ZSR relationship is used to set water levels. Units are
         !  10**(-6) kg water/ (cubic meter of air)
         !  start with ammomium sulfate solution without nitrate
       
         CALL AWATER_new( IRH, TSO4, YNH4, TNO3, AH2O ) !**** note TNO3
         WH2O = 1.0d-3 * AH2O
         ASO4 = TSO4   * MWSO4

         ! In sulfate poor case, Sulfate ion is preferred
         ! Set bisulfate equal to zero [rjp, 12/12/01]
         AHSO4 = 0.0d0
         ANO3  = 0.0d0
         ANH4  = YNH4 * MWNH4
         WFRAC = AH2O / ( ASO4 + ANH4 +  AH2O )

        !IF ( WFRAC .EQ. 0.0 )  RETURN   ! No water
        IF ( WFRAC .LT. 0.2d0 ) THEN

           ! "dry" ammonium sulfate and ammonium nitrate
           ! compute free ammonia 
           FNH3 = TNH4 - TWOSO4
           CC   = TNO3 * FNH3 - K3

           ! check for not enough to support aerosol
           IF ( CC .LE. 0.0d0 ) THEN
              XNO3 = 0.0d0
           ELSE
              AA   = 1.0d0
              BB   = -( TNO3 + FNH3 )
              DISC = BB * BB - 4.0d0 * CC

              ! Check for complex roots of the quadratic
              ! set retain initial values of nitrate and RETURN 
              ! if complex roots are found
              IF ( DISC .LT. 0.0d0 ) THEN
                 XNO3  = 0.0d0
                 AH2O  = 1000.0d0 * WH2O
                 YNH4  = TWOSO4
                 ASO4  = TSO4 * MWSO4
                 AHSO4 = 0.0d0
                 ANH4  = YNH4 * MWNH4
                 GNH3  = MWNH3 * MAX( FLOOR, ( TNH4 - YNH4 ) )
                 GNO3  = GNO3_IN
                 ANO3  = ANO3_IN
                 RETURN
              ENDIF

              ! to get here, BB .lt. 0.0, CC .gt. 0.0 always
              DD  = SQRT( DISC )
              XXQ = -0.5d0 * ( BB + SIGN ( 1.0, BB ) * DD )


              ! Since both roots are positive, select smaller root.
              XNO3 = MIN( XXQ / AA, CC / XXQ )

           ENDIF                ! CC .LE. 0.0

           AH2O  = 1000.0d0 * WH2O
           YNH4  = TWOSO4 + XNO3          
           ASO4  = TSO4 * MWSO4
           AHSO4 = FLOOR
           ANO3  = XNO3 * MWNO3
           ANH4  = YNH4 * MWNH4
           GNH3  = MWNH3 * MAX( FLOOR, ( TNH4 - YNH4 )  )      
           GNO3  = MAX( FLOOR, ( TMASSHNO3 - ANO3 ) )    
           RETURN
        ENDIF                  ! WFRAC .LT. 0.2

        ! liquid phase containing completely neutralized sulfate and
        ! some nitrate.  Solve for composition and quantity.
        MAS    = TSO4 / WH2O
        MAN    = 0.0d0
        XNO3   = 0.0d0
        YNH4   = TWOSO4
        PHIOLD = 1.0d0

        !===============================================================
        ! Start loop for iteration
        !
        ! The assumption here is that all sulfate is ammonium sulfate,
        ! and is supersaturated at lower relative humidities.
        !===============================================================
        DO NNN = 1, 50 ! loop count reduced 0409/2001 by FSB

           NITR  = NNN
           GASQD = GAMAAN * GAMAAN
           WSQD  = WH2O * WH2O
           KW2   = KAN * WSQD / GASQD
           AA    = 1.0 - KW2
           BB    = TWOSO4 + KW2 * ( TNO3 + TNH4 - TWOSO4 )
           CC    = -KW2 * TNO3 * ( TNH4 - TWOSO4 )

           ! This is a quadratic for XNO3 [MICROMOLES / M**3] 
           ! of nitrate in solution
           DISC = BB * BB - 4.0d0 * AA * CC

           ! Check for complex roots, retain inital values and RETURN
           IF ( DISC .LT. 0.0 ) THEN
              XNO3  = 0.0d0
              AH2O  = 1000.0d0 * WH2O
              YNH4  = TWOSO4
              ASO4  = TSO4 * MWSO4
              AHSO4 = FLOOR     ! [rjp, 12/12/01]
              ANH4  = YNH4 * MWNH4
              GNH3  = MWNH3 * MAX( FLOOR, (TNH4 - YNH4 ) )
              GNO3  = GNO3_IN
              ANO3  = ANO3_IN
              RETURN
           ENDIF
           
           ! Deal with degenerate case (yoj)
           IF ( AA .NE. 0.0d0 ) THEN
              DD  = SQRT( DISC )
              XXQ = -0.5d0 * ( BB + SIGN( 1.0, BB ) * DD )
              RR1 = XXQ / AA
              RR2 = CC / XXQ

              ! choose minimum positve root
              IF ( ( RR1 * RR2 ) .LT. 0.0d0 ) THEN
                 XNO3 = MAX( RR1, RR2 )
              ELSE
                 XNO3 = MIN( RR1, RR2 )
              ENDIF
           ELSE
              XNO3 = - CC / BB  ! AA equals zero here.
           ENDIF

           XNO3 = MIN( XNO3, TNO3 )
           
           ! This version assumes no solid sulfate forms (supersaturated )
           ! Now update water
           CALL AWATER_new ( IRH, TSO4, YNH4, XNO3, AH2O )

           ! ZSR relationship is used to set water levels. Units are
           ! 10**(-6) kg water/ (cubic meter of air).  The conversion 
           ! from micromoles to moles is done by the units of WH2O.
           WH2O = 1.0d-3 * AH2O 

           ! Ionic balance determines the ammonium in solution.
           MAN  = XNO3 / WH2O
           MAS  = TSO4 / WH2O
           MNH4 = 2.0d0 * MAS + MAN
           YNH4 = MNH4 * WH2O

           ! MAS, MAN and MNH4 are the aqueous concentrations of sulfate, 
           ! nitrate, and ammonium in molal units (moles/(kg water) ).
           STION    = 3.0d0 * MAS + MAN
           CAT( 1 ) = 0.0d0
           CAT( 2 ) = MNH4
           AN ( 1 ) = MAS
           AN ( 2 ) = MAN
           AN ( 3 ) = 0.0d0
           CALL ACTCOF_new ( CAT, AN, GAMS, MOLNU, PHIBAR )
           GAMAAN = GAMS( 2, 2 )

           ! Use GAMAAN for convergence control
           EROR   = ABS( GAMOLD - GAMAAN ) / GAMOLD
           GAMOLD = GAMAAN

           ! Check to see if we have a solution
           IF ( EROR .LE. TOLER1 ) THEN
              ASO4  = TSO4 * MWSO4
              AHSO4 = 0.0d0       ! [rjp, 12/12/01]
              ANO3  = XNO3 * MWNO3
              ANH4  = YNH4 * MWNH4
              GNO3  = MAX( FLOOR, ( TMASSHNO3  - ANO3 ) )
              GNH3  = MWNH3 * MAX( FLOOR, ( TNH4 - YNH4 ) )
              AH2O  = 1000.0d0 * WH2O
              RETURN
           ENDIF

        ENDDO

        ! If after NITR iterations no solution is found, then:
        ! FSB retain the initial values of nitrate particle and vapor
        ! note whether or not convert all bisulfate to sulfate
        ASO4  = TSO4 * MWSO4
        AHSO4 = FLOOR      
        XNO3  = TNO3 / MWNO3
        YNH4  = TWOSO4
        ANH4  = YNH4 * MWNH4

        CALL AWATER_new ( IRH, TSO4, YNH4, XNO3, AH2O )

        GNO3  = GNO3_IN
        ANO3  = ANO3_IN
        GNH3  = MAX( FLOOR, MWNH3 * (TNH4 - YNH4 ) )
        RETURN

      !================================================================
      ! Low Ammonia Case 
      !
      ! Coded by Dr. Francis S. Binkowski 12/8/91.(4/26/95)
      ! modified 8/28/98
      ! modified 04/09/2001
      !       
      ! All cases covered by this logic
      !=================================================================
      ELSE

         WH2O = 0.0d0
         CALL AWATER_new ( IRH, TSO4, TNH4, TNO3, AH2O )
         WH2O = 1.0d-3 * AH2O
         ZH2O = AH2O

         ! convert 10**(-6) kg water/(cubic meter of air) to micrograms 
         ! of water per cubic meter of air (1000 g = 1 kg)
         ! in sulfate rich case, preferred form is HSO4-
         !ASO4 = TSO4 * MWSO4
         ASO4  = FLOOR          ![rjp, 12/12/01]
         AHSO4 = TSO4 * MWSO4   ![rjp, 12/12/01]
         ANH4  = TNH4 * MWNH4
         ANO3  = ANO3_IN
         GNO3  = TMASSHNO3 - ANO3
         GNH3  = FLOOR 
        
         !==============================================================
         ! *** Examine special cases and return if necessary.
         !         
         ! FSB For values of RATIO less than 0.5 do no further 
         ! calculations.  The code will cycle and still predict the 
         ! same amount of ASO4, ANH4, ANO3, AH2O so terminate early 
         ! to swame computation
         !==============================================================
         IF ( RATIO .LT. 0.5d0 ) RETURN ! FSB 04/09/2001 
    
         ! Check for zero water.
         IF ( WH2O .EQ. 0.0d0 ) RETURN
         ZSO4 = TSO4 / WH2O

         ! ZSO4 is the molality of total sulfate i.e. MSO4 + MHSO4
         ! do not solve for aerosol nitrate for total sulfate molality
         ! greater than 11.0 because the model parameters break down
         !### IF ( ZSO4 .GT. 11.0 ) THEN
         IF ( ZSO4 .GT. 9.0 ) THEN ! 18 June 97
            RETURN
         ENDIF

         ! *** Calculation may now proceed.
         !
         ! First solve with activity coeffs of 1.0, then iterate.
         PHIOLD = 1.0d0
         GAMANA = 1.0d0
         GAMAS1 = 1.0d0
         GAMAS2 = 1.0d0
         GAMAAB = 1.0d0
         GAMOLD = 1.0d0

         ! All ammonia is considered to be aerosol ammonium.
         MNH4 = TNH4 / WH2O

         ! MNH4 is the molality of ammonium ion.
         YNH4 = TNH4

         ! loop for iteration
         DO NNN = 1, 50    ! loop count reduced 04/09/2001 by FSB
            NITR = NNN

            !------------------------------------------------------------
            ! Add robustness: now check if GAMANA or GAMAS1 is too small
            ! for the division in RKNA and RK2SA. If they are, return w/ 
            ! original values: basically replicate the procedure used 
            ! after the current DO-loop in case of no-convergence
            ! (phs, bmy, rjp, 4/10/08)
            ! Now uses IS_SAFE_DIV to avoid compiler/machine dependency 
            ! and to check for both underlow and overflow. Also 
            ! use REAL4 flag to avoid under/overflow when computing A0 
            ! and A1 from RKNA and RK2SA (phs, 5/28/08)
            !------------------------------------------------------------
            IF ( .NOT. (&
                IS_SAFE_DIV( GAMAS2, GAMAS1*GAMAS1*GAMAS1, R4=.TRUE. )& 
                .AND. IS_SAFE_DIV( KNA, GAMANA*GAMANA, R4=.TRUE. ) &
               ) ) THEN
               
               WRITE(6,*) 'RPMARES: not safe to divide...exit'
               CALL flush(6)
               GOTO 100
            ENDIF
            
            ! set up equilibrium constants including activities
            ! solve the system for hplus first then sulfate & nitrate
            RK2SA  = K2SA * GAMAS2 * GAMAS2 / (GAMAS1 * GAMAS1 * GAMAS1)
            RKNA   = KNA / ( GAMANA * GAMANA )
            RKNWET = RKNA * WH2O
            T21    = ZSO4 - MNH4
            T221   = ZSO4 + T21

            ! set up coefficients for cubic
            A2 = RK2SA + RKNWET - T21
            A1 = RK2SA * RKNWET - T21 * ( RK2SA + RKNWET )&
                - RK2SA * ZSO4 - RKNA * TNO3
            A0 = - (T21 * RK2SA * RKNWET&
                + RK2SA * RKNWET * ZSO4 + RK2SA * RKNA * TNO3 )

            CALL CUBIC ( A2, A1, A0, NR, CRUTES )

            ! Code assumes the smallest positive root is in CRUTES(1)
            ! But, it can be negative (see CUBIC, case of one real root, 
            ! but can also be propagated by over/underflown)... if it is
            ! the case then return with original values (phs, 5/27/08)
            HPLUS = CRUTES( 1 )
            IF (HPLUS <= 0d0) GOTO 100
            BAL   = HPLUS **3 + A2 * HPLUS**2 + A1 * HPLUS + A0

            ! molality of sulfate ion
            MSO4  = RK2SA * ZSO4 / ( HPLUS + RK2SA ) 

            ! molality of bisulfate ion
            ! MAX added 04/09/2001 by FSB
            MHSO4 = MAX( 1.0d-10, ZSO4 - MSO4 ) 

            ! molality of nitrate ion
            MNA   = RKNA * TNO3 / ( HPLUS + RKNWET ) 
            MNA   = MAX( 0.0d0, MNA )
            MNA   = MIN( MNA, TNO3 / WH2O )
            XNO3  = MNA * WH2O
            ANO3  = MNA * WH2O * MWNO3
            GNO3  = MAX( FLOOR, TMASSHNO3 - ANO3 )
            ASO4  = MSO4 * WH2O * MWSO4 ![rjp, 12/12/01]
            AHSO4 = MHSO4 * WH2O * MWSO4 ![rjp, 12/12/01]

            ! Calculate ionic strength
            STION = 0.5d0 * ( HPLUS + MNA + MNH4 + MHSO4 + 4.0d0 * MSO4)

            ! Update water
            CALL AWATER_new ( IRH, TSO4, YNH4, XNO3, AH2O )

            ! Convert 10**(-6) kg water/(cubic meter of air) to micrograms 
            ! of water per cubic meter of air (1000 g = 1 kg)
            WH2O     = 1.0d-3 * AH2O
            CAT( 1 ) = HPLUS
            CAT( 2 ) = MNH4
            AN ( 1 ) = MSO4
            AN ( 2 ) = MNA
            AN ( 3 ) = MHSO4

            CALL ACTCOF_new ( CAT, AN, GAMS, MOLNU, PHIBAR )

            GAMANA = GAMS( 1, 2 )
            GAMAS1 = GAMS( 1, 1 )
            GAMAS2 = GAMS( 1, 3 )
            GAMAAN = GAMS( 2, 2 )

            ! NOTE: Improved for robustness!
            GAMAHAT = ( GAMAS2 * GAMAS2 / ( GAMAAB * GAMAAB ) )
            BHAT = KHAT * GAMAHAT
            !### EROR = ABS ( ( PHIOLD - PHIBAR ) / PHIOLD )
            !### PHIOLD = PHIBAR
            EROR = ABS ( GAMOLD - GAMAHAT ) / GAMOLD
            GAMOLD = GAMAHAT

            ! return with good solution
            IF ( EROR .LE. TOLER2 ) THEN
               RETURN
            ENDIF
            
         ENDDO

         ! after NITR iterations, failure to solve the system
         ! convert all ammonia to aerosol ammonium and return input
         ! values of NO3 and HNO3
 100     ANH4 = TNH4 * MWNH4
         GNH3 = FLOOR
         GNO3 = GNO3_IN
         ANO3 = ANO3_IN
         ASO4 = TSO4 * MWSO4    ! [rjp, 12/17/01]
         AHSO4= FLOOR           ! [rjp, 12/17/01]

         CALL AWATER_new ( IRH, TSO4, TNH4, TNO3, AH2O )
         
         RETURN

      ENDIF                     ! ratio .gt. 2.0

      ! Return to calling program
      END SUBROUTINE RPMARES_new

      SUBROUTINE AWATER_new( IRHX, MSO4, MNH4, MNO3, WH2O )
!
!******************************************************************************
! NOTE!!! wh2o is returned in micrograms / cubic meter
!         mso4,mnh4,mno3 are in microMOLES / cubic meter
!
!  This  version uses polynomials rather than tables, and uses empirical
! polynomials for the mass fraction of solute (mfs) as a function of water
! activity
!   where:
!
!            mfs = ms / ( ms + mw)
!             ms is the mass of solute
!             mw is the mass of water.
!
!  Define y = mw/ ms
!
!  then  mfs = 1 / (1 + y)
!
!    y can then be obtained from the values of mfs as
!
!             y = (1 - mfs) / mfs
!
!
!     the aerosol is assumed to be in a metastable state if the rh is
!     is below the rh of deliquescence, but above the rh of crystallization.
!
!     ZSR interpolation is used for sulfates with x ( the molar ratio of
!     ammonium to sulfate in eh range 0 <= x <= 2, by sections.
!     section 1: 0 <= x < 1
!     section 2: 1 <= x < 1.5
!     section 3: 1.5 <= x < 2.0
!     section 4: 2 <= x
!     In sections 1 through 3, only the sulfates can affect the amount of water
!     on the particles.
!     In section 4, we have fully neutralized sulfate, and extra ammonium which
!     allows more nitrate to be present. Thus, the ammount of water is
!     calculated
!     using ZSR for ammonium sulfate and ammonium nitrate. Crystallization is
!     assumed to occur in sections 2,3,and 4. See detailed discussion below.
!
!
!
! definitions:
!     mso4, mnh4, and mno3 are the number of micromoles/(cubic meter of air)
!      for sulfate, ammonium, and nitrate respectively
!     irhx is the relative humidity (%)
!     wh2o is the returned water amount in micrograms / cubic meter of air
!     x is the molar ratio of ammonium to sulfate
!     y0,y1,y1.5, y2 are the water contents in mass of water/mass of solute
!     for pure aqueous solutions with x equal 1, 1.5, and 2 respectively.
!     y3 is the value of the mass ratio of water to solute for
!     a pure ammonium nitrate  solution.
!
!
!     coded by Dr. Francis S. Binkowski, 4/8/96.
!
! *** modified 05/30/2000 
!     The use of two values of mfs at an ammonium to sulfate ratio 
!     representative of ammonium sulfate led to an minor inconsistancy 
!     in nitrate behavior as the ratio went from a value less than two
!     to a value greater than two and vice versa with either ammonium 
!     held constant and sulfate changing, or sulfate held constant and 
!     ammonium changing. the value of Chan et al. (1992) is the only value
!     now used. 
!
! *** Modified 09/25/2002
!     Ported into "rpmares_mod.f".  Now declare all variables with REAL*8.
!     Also cleaned up comments and made cosmetic changes.  Force double 
!     precision explicitly with "D" exponents. 
!******************************************************************************
!
      ! Arguments
      INTEGER           :: IRHX
      REAL*8            :: MSO4, MNH4, MNO3, WH2O

      ! Local variables
      INTEGER           :: IRH
      REAL*8            :: TSO4,  TNH4,  TNO3,  X,      AW,     AWC
      REAL*8            :: MFS0,  MFS1,  MFS15, Y , MFS2
      REAL*8            :: Y0,    Y1,    Y15,   Y2,     Y3,     Y40 
      REAL*8            :: Y140,  Y1540, YC,    MFSSO4, MFSNO3

      ! Molecular weight parameters
      REAL*8, PARAMETER :: MWSO4  = 96.0636d0
      REAL*8, PARAMETER :: MWNH4  = 18.0985d0
      REAL*8, PARAMETER :: MWNO3  = 62.0649d0
      REAL*8, PARAMETER :: MW2    = MWSO4 + 2.0d0 * MWNH4
      REAL*8, PARAMETER :: MWANO3 = MWNO3 + MWNH4 
      
      !=================================================================
      ! The polynomials use data for aw as a function of mfs from Tang 
      ! and Munkelwitz, JGR 99: 18801-18808, 1994.  The polynomials were 
      ! fit to Tang's values of water activity as a function of mfs.
      ! 
      ! *** coefficients of polynomials fit to Tang and Munkelwitz data
      !     now give mfs as a function of water activity.
      !=================================================================
      REAL*8 :: C1(4)  = (/ 0.9995178d0,  -0.7952896d0, &
                           0.99683673d0, -1.143874d0 /)

      REAL*8 :: C15(4) = (/ 1.697092d0, -4.045936d0, &
                           5.833688d0, -3.463783d0 /)
      
      REAL*8 :: C2(4)  = (/ 2.085067d0, -6.024139d0, &
                           8.967967d0, -5.002934d0 /)

      !=================================================================
      ! The following coefficients are a fit to the data in Table 1 of
      !    Nair & Vohra, J. Aerosol Sci., 6: 265-271, 1975
      !      data c0/0.8258941, -1.899205, 3.296905, -2.214749 /
      !
      ! New data fit to data from
      !       Nair and Vohra J. Aerosol Sci., 6: 265-271, 1975
      !       Giaque et al. J.Am. Chem. Soc., 82: 62-70, 1960
      !       Zeleznik J. Phys. Chem. Ref. Data, 20: 157-1200
      !=================================================================
      REAL*8 :: C0(4)  =  (/ 0.798079d0, -1.574367d0, &
                            2.536686d0, -1.735297d0 /)

      !=================================================================
      ! Polynomials for ammonium nitrate and ammonium sulfate are from:
      ! Chan et al.1992, Atmospheric Environment (26A): 1661-1673.
      !=================================================================
      REAL*8 :: KNO3(6) = (/  0.2906d0,   6.83665d0, -26.9093d0,&
                            46.6983d0, -38.803d0,    11.8837d0 /)
      
      REAL*8 :: KSO4(6) = (/   2.27515d0, -11.147d0,   36.3369d0,&
                            -64.2134d0,   56.8341d0, -20.0953d0 /)

      !=================================================================
      ! AWATER begins here!
      !=================================================================
      
      ! Check range of per cent relative humidity
      IRH  = IRHX
      IRH  = MAX( 1, IRH )
      IRH  = MIN( IRH, 100 )

      ! Water activity = fractional relative humidity
      AW   = DBLE( IRH ) / 100.0d0
      TSO4 = MAX( MSO4 , 0.0d0 )
      TNH4 = MAX( MNH4 , 0.0d0 )
      TNO3 = MAX( MNO3 , 0.0d0 )
      X    = 0.0d0

      ! If there is non-zero sulfate calculate the molar ratio
      ! otherwise check for non-zero nitrate and ammonium
      IF ( TSO4 .GT. 0.0d0 ) THEN
         X = TNH4 / TSO4
      ELSE
         IF ( TNO3 .GT. 0.0d0 .AND. TNH4 .GT. 0.0d0 ) X = 10.0d0
      ENDIF

      ! *** begin screen on x for calculating wh2o
      IF ( X .LT. 1.0d0 ) THEN
         MFS0 = POLY4( C0, AW )
         MFS1 = POLY4( C1, AW )
         Y0   = ( 1.0d0 - MFS0 ) / MFS0
         Y1   = ( 1.0d0 - MFS1 ) / MFS1
         Y    = ( 1.0d0 - X    ) * Y0 + X * Y1

      ELSE IF ( X .LT. 1.5d0 ) THEN

         IF ( IRH .GE. 40 ) THEN
            MFS1  = POLY4( C1,  AW )
            MFS15 = POLY4( C15, AW )
            Y1    = ( 1.0d0 - MFS1  ) / MFS1
            Y15   = ( 1.0d0 - MFS15 ) / MFS15
            Y     = 2.0d0 * ( Y1 * ( 1.5d0 - X ) + Y15 *( X - 1.0d0 ) )

         !==============================================================
         ! Set up for crystalization
         !
         ! Crystallization is done as follows:
         !
         ! For 1.5 <= x, crystallization is assumed to occur 
         ! at rh = 0.4
         !
         ! For x <= 1.0, crystallization is assumed to occur at an 
         ! rh < 0.01, and since the code does not allow ar rh < 0.01, 
         ! crystallization is assumed not to occur in this range.
         !
         ! For 1.0 <= x <= 1.5 the crystallization curve is a straignt 
         ! line from a value of y15 at rh = 0.4 to a value of zero at 
         ! y1. From point B to point A in the diagram.  The algorithm 
         ! does a double interpolation to calculate the amount of
         ! water.
         !
         !        y1(0.40)               y15(0.40)
         !         +                     + Point B
         !
         !
         !
         !
         !         +--------------------+
         !       x=1                   x=1.5
         !      Point A
         !==============================================================
         ELSE

            ! rh along the crystallization curve.
            AWC = 0.80d0 * ( X - 1.0d0 ) 
            Y   = 0.0d0

            ! interpolate using crystalization curve
            IF ( AW .GE. AWC ) THEN 
               MFS1  = POLY4( C1,  0.40 )
               MFS15 = POLY4( C15, 0.40 )
               Y140  = ( 1.0d0 - MFS1  ) / MFS1
               Y1540 = ( 1.0d0 - MFS15 ) / MFS15
               Y40   = 2.0d0 * ( Y140  * ( 1.5d0 - X ) + &
                                Y1540 * ( X - 1.0d0 ) )

               ! Y along crystallization curve
               YC   = 2.0d0 * Y1540 * ( X - 1.0d0 ) 
               Y    = Y40 - (Y40 - YC) * (0.40d0 - AW) / (0.40d0 - AWC)
            ENDIF              
         ENDIF                 

      ELSE IF ( X .LT. 2.0d0 ) then               ! changed 12/11/2000 by FSB
         Y = 0.0D0

         IF ( IRH .GE. 40 ) THEN 
            MFS15  = POLY4( C15, AW )            
!            MFS2  = POLY4( C2,  AW )
            Y15    = ( 1.0d0 - MFS15 ) / MFS15
!            y2    = ( 1.0d0 - MFS2  ) / MFS2
            MFSSO4 = POLY6( KSO4, AW )             ! Changed 05/30/2000 by FSB
            Y2     = ( 1.0d0 - MFSSO4 ) / MFSSO4
            Y      = 2.0d0 * (Y15 * (2.0d0 - X) + Y2 * (X - 1.5d0) )
         ENDIF                  
      ELSE                                 ! 2.0 <= x changed 12/11/2000 by FSB

         !==============================================================
         ! Regime where ammonium sulfate and ammonium nitrate are 
         ! in solution.
         ! 
         ! following cf&s for both ammonium sulfate and ammonium nitrate
         ! check for crystallization here. their data indicate a 40% 
         ! value is appropriate.
         !==============================================================
         Y2 = 0.0d0
         Y3 = 0.0d0

         IF ( IRH .GE. 40 ) THEN
            MFSSO4 = POLY6( KSO4, AW )
            MFSNO3 = POLY6( KNO3, AW )
            Y2     = ( 1.0d0 - MFSSO4 ) / MFSSO4
            Y3     = ( 1.0d0 - MFSNO3 ) / MFSNO3

         ENDIF

      ENDIF                     ! end of checking on x

      !=================================================================
      ! Now set up output of WH2O
      ! WH2O units are micrograms (liquid water) / cubic meter of air
      !=================================================================
      IF ( X .LT. 2.0D0 ) THEN  ! changed 12/11/2000 by FSB

         WH2O =  Y * ( TSO4 * MWSO4 + MWNH4 * TNH4 )

      ELSE

         ! this is the case that all the sulfate is ammonium sulfate
         ! and the excess ammonium forms ammonum nitrate
         WH2O =   Y2 * TSO4 * MW2 + Y3 * TNO3 * MWANO3

      ENDIF

      ! Return to calling program
      END SUBROUTINE AWATER_new


       SUBROUTINE ACTCOF_new( CAT, AN, GAMA, MOLNU, PHIMULT )
!
!******************************************************************************
!
! DESCRIPTION:
!
!  This subroutine computes the activity coefficients of (2NH4+,SO4--),
!  (NH4+,NO3-),(2H+,SO4--),(H+,NO3-),AND (H+,HSO4-) in aqueous
!  multicomponent solution, using Bromley's model and Pitzer's method.
!
! REFERENCES:
!
!   Bromley, L.A. (1973) Thermodynamic properties of strong electrolytes
!     in aqueous solutions.  AIChE J. 19, 313-320.
!
!   Chan, C.K. R.C. Flagen, & J.H.  Seinfeld (1992) Water Activities of
!     NH4NO3 / (NH4)2SO4 solutions, Atmos. Environ. (26A): 1661-1673.
!
!   Clegg, S.L. & P. Brimblecombe (1988) Equilibrium partial pressures
!     of strong acids over saline solutions - I HNO3,
!     Atmos. Environ. (22): 91-100
!
!   Clegg, S.L. & P. Brimblecombe (1990) Equilibrium partial pressures
!     and mean activity and osmotic coefficients of 0-100% nitric acid
!     as a function of temperature,   J. Phys. Chem (94): 5369 - 5380
!
!   Pilinis, C. and J.H. Seinfeld (1987) Continued development of a
!     general equilibrium model for inorganic multicomponent atmospheric
!     aerosols.  Atmos. Environ. 21(11), 2453-2466.
!
!
!
!
! ARGUMENT DESCRIPTION:
!
!     CAT(1) : conc. of H+    (moles/kg)
!     CAT(2) : conc. of NH4+  (moles/kg)
!     AN(1)  : conc. of SO4-- (moles/kg)
!     AN(2)  : conc. of NO3-  (moles/kg)
!     AN(3)  : conc. of HSO4- (moles/kg)
!     GAMA(2,1)    : mean molal ionic activity coeff for (2NH4+,SO4--)
!     GAMA(2,2)    :  "    "     "       "       "    "  (NH4+,NO3-)
!     GAMA(2,3)    :  "    "     "       "       "    "  (NH4+. HSO4-)
!     GAMA(1,1)    :  "    "     "       "       "    "  (2H+,SO4--)
!     GAMA(1,2)    :  "    "     "       "       "    "  (H+,NO3-)
!     GAMA(1,3)    :  "    "     "       "       "    "  (H+,HSO4-)
!     MOLNU   : the total number of moles of all ions.
!     PHIMULT : the multicomponent paractical osmotic coefficient.
!
! REVISION HISTORY:
!      Who       When        Detailed description of changes
!   ---------   --------  -------------------------------------------
!   S.Roselle   7/26/89   Copied parts of routine BROMLY, and began this
!                         new routine using a method described by Pilinis
!                         and Seinfeld 1987, Atmos. Envirn. 21 pp2453-2466.
!   S.Roselle   7/30/97   Modified for use in Models-3
!   F.Binkowski 8/7/97    Modified coefficients BETA0, BETA1, CGAMA
!   R.Yantosca  9/25/02   Ported into "rpmares_mod.f" for GEOS-CHEM.  Cleaned
!                         up comments, etc.  Also force double precision by
!                         declaring REALs as REAL*8 and by using "D" exponents.
!******************************************************************************
!

      ! Error codes


     
      !=================================================================
      ! PARAMETERS and their descriptions:
      !=================================================================
      INTEGER, PARAMETER :: NCAT = 2         ! number of cation
      INTEGER, PARAMETER :: NAN  = 3         ! number of anions
      REAL*8,  PARAMETER :: XSTAT0 = 0       ! Normal, successful completion
      REAL*8,  PARAMETER :: XSTAT1 = 1       ! File I/O error
      REAL*8,  PARAMETER :: XSTAT2 = 2       ! Execution error
      REAL*8,  PARAMETER :: XSTAT3 = 3       ! Special  error

      !=================================================================
      ! ARGUMENTS and their descriptions
      !=================================================================
      REAL*8             :: MOLNU            ! tot # moles of all ions
      REAL*8             :: PHIMULT          ! multicomponent paractical 
                                             !   osmotic coef
      REAL*8             :: CAT(NCAT)        ! cation conc in moles/kg (input)
      REAL*8             :: AN(NAN)          ! anion conc in moles/kg (input)
      REAL*8             :: GAMA(NCAT,NAN)   ! mean molal ionic activity coefs

      !=================================================================
      ! SCRATCH LOCAL VARIABLES and their descriptions:
      !=================================================================
      INTEGER            :: IAN              ! anion indX
      INTEGER            :: ICAT             ! cation indX
      REAL*8             :: FGAMA            !
      REAL*8             :: I                ! ionic strength
      REAL*8             :: R                !
      REAL*8             :: S                !
      REAL*8             :: TA               !
      REAL*8             :: TB               !
      REAL*8             :: TC               !
      REAL*8             :: TEXPV            !
      REAL*8             :: TRM              !
      REAL*8             :: TWOI             ! 2*ionic strength
      REAL*8             :: TWOSRI           ! 2*sqrt of ionic strength
      REAL*8             :: ZBAR             !
      REAL*8             :: ZBAR2            !
      REAL*8             :: ZOT1             !
      REAL*8             :: SRI              ! square root of ionic strength
      REAL*8             :: F2(NCAT)         !
      REAL*8             :: F1(NAN)          !
      REAL*8             :: BGAMA (NCAT,NAN) !
      REAL*8             :: X     (NCAT,NAN) !
      REAL*8             :: M     (NCAT,NAN) ! molality of each electrolyte
      REAL*8             :: LGAMA0(NCAT,NAN) ! binary activity coefficients
      REAL*8             :: Y     (NAN,NCAT) !
      REAL*8             :: BETA0 (NCAT,NAN) ! binary activity coef parameter
      REAL*8             :: BETA1 (NCAT,NAN) ! binary activity coef parameter
      REAL*8             :: CGAMA (NCAT,NAN) ! binary activity coef parameter
      REAL*8             :: V1    (NCAT,NAN) ! # of cations in electrolyte
                                             !   formula
      REAL*8             :: V2    (NCAT,NAN) ! # of anions in electrolyte
                                             !   formula
      ! absolute value of charges of cation
      REAL*8             :: ZP(NCAT) = (/ 1.0d0, 1.0d0 /)         

      ! absolute value of charges of anion
      REAL*8             :: ZM(NAN)  = (/ 2.0d0, 1.0d0, 1.0d0 /)         

      ! Character values.
      CHARACTER(LEN=120)      :: XMSG  = ' '
      CHARACTER(LEN=16), SAVE :: PNAME = ' driver program name'

      !================================================================
      ! *** Sources for the coefficients BETA0, BETA1, CGAMA
      ! (1,1);(1,3)  - Clegg & Brimblecombe (1988)
      ! (2,3)        - Pilinis & Seinfeld (1987), cgama different
      ! (1,2)        - Clegg & Brimblecombe (1990)
      ! (2,1);(2,2)  - Chan, Flagen & Seinfeld (1992)
      !================================================================
 
      ! now set the basic constants, BETA0, BETA1, CGAMA
      DATA BETA0(1,1) /2.98d-2/,      BETA1(1,1) / 0.0d0/,&
          CGAMA(1,1) /4.38d-2/                                 ! 2H+SO4-

      DATA BETA0(1,2) /  1.2556d-1/,  BETA1(1,2) / 2.8778d-1/,&
          CGAMA(1,2) / -5.59d-3/                               ! HNO3

      DATA BETA0(1,3) / 2.0651d-1/,   BETA1(1,3) / 5.556d-1/,&
          CGAMA(1,3) /0.0d0/                                   ! H+HSO4-

      DATA BETA0(2,1) / 4.6465d-2/,   BETA1(2,1) /-0.54196d0/,&
          CGAMA(2,1) /-1.2683d-3/                              ! (NH4)2SO4

      DATA BETA0(2,2) /-7.26224d-3/,  BETA1(2,2) /-1.168858d0/,&
          CGAMA(2,2) / 3.51217d-5/                             ! NH4NO3

      DATA BETA0(2,3) / 4.494d-2/,    BETA1(2,3) / 2.3594d-1/,&
          CGAMA(2,3) /-2.962d-3/                               ! NH4HSO4

      DATA V1(1,1), V2(1,1) / 2.0d0, 1.0d0 /     ! 2H+SO4-
      DATA V1(2,1), V2(2,1) / 2.0d0, 1.0d0 /     ! (NH4)2SO4
      DATA V1(1,2), V2(1,2) / 1.0d0, 1.0d0 /     ! HNO3
      DATA V1(2,2), V2(2,2) / 1.0d0, 1.0d0 /     ! NH4NO3
      DATA V1(1,3), V2(1,3) / 1.0d0, 1.0d0 /     ! H+HSO4-
      DATA V1(2,3), V2(2,3) / 1.0d0, 1.0d0 /     ! NH4HSO4

      !=================================================================
      ! ACTCOF begins here!
      !=================================================================

      ! Compute ionic strength
      I = 0.0d0

      DO ICAT = 1, NCAT
         I = I + CAT( ICAT ) * ZP( ICAT ) * ZP( ICAT )
      ENDDO

      DO IAN = 1, NAN
         I = I + AN( IAN ) * ZM( IAN ) * ZM( IAN )
      ENDDO

      I = 0.5d0 * I

      ! check for problems in the ionic strength
      IF ( I .EQ. 0.0d0 ) THEN

         DO IAN  = 1, NAN
         DO ICAT = 1, NCAT
            GAMA( ICAT, IAN ) = 0.0d0
         ENDDO
         ENDDO

         XMSG = 'Ionic strength is zero...returning zero activities'
         !CALL M3WARN ( PNAME, 0, 0, XMSG )
         RETURN

      ELSE IF ( I .LT. 0.0d0 ) THEN
         XMSG = 'Ionic strength below zero...negative concentrations'
         write(6,*)xmsg
         call flush(6)
         !CALL M3EXIT ( PNAME, 0, 0, XMSG, XSTAT1 )
      ENDIF

      ! Compute some essential expressions
      SRI    = SQRT( I )
      TWOSRI = 2.0d0 * SRI
      TWOI   = 2.0d0 * I
      TEXPV  = 1.0d0 - EXP( -TWOSRI ) * ( 1.0d0 + TWOSRI - TWOI )
      R      = 1.0d0 + 0.75d0 * I
      S      = 1.0d0 + 1.5d0  * I
      ZOT1   = 0.511d0 * SRI / ( 1.0d0 + SRI )

      ! Compute binary activity coeffs
      FGAMA = -0.392d0 * ( ( SRI / ( 1.0d0 + 1.2d0 * SRI )&
           + ( 2.0d0 / 1.2d0 ) * LOG( 1.0d0 + 1.2d0 * SRI ) ) )

      DO ICAT = 1, NCAT
      DO IAN  = 1, NAN

         BGAMA( ICAT, IAN ) = 2.0d0 * BETA0( ICAT, IAN )&
             + ( 2.0d0 * BETA1( ICAT, IAN ) / ( 4.0d0 * I ) )&
             * TEXPV

         ! Compute the molality of each electrolyte for given ionic strength
         M( ICAT, IAN ) = ( CAT( ICAT )**V1( ICAT, IAN )&
                        *   AN( IAN )**V2( ICAT, IAN ) )**( 1.0d0&
                        / ( V1( ICAT, IAN ) + V2( ICAT, IAN ) ) )

         ! Calculate the binary activity coefficients
         LGAMA0( ICAT, IAN ) = ( ZP( ICAT ) * ZM( IAN ) * FGAMA&
             + M( ICAT, IAN )&
             * ( 2.0d0 * V1( ICAT, IAN ) * V2( ICAT, IAN )&
             / ( V1( ICAT, IAN ) + V2( ICAT, IAN ) )&
             * BGAMA( ICAT, IAN ) )                 &
             + M( ICAT, IAN ) * M( ICAT, IAN )      &
             * ( 2.0d0 * ( V1( ICAT, IAN )          &
             * V2( ICAT, IAN ) )**1.5d0             &
             / ( V1( ICAT, IAN ) + V2( ICAT, IAN ) )&
             * CGAMA( ICAT, IAN ) ) ) / 2.302585093d0

      ENDDO
      ENDDO

      ! prepare variables for computing the multicomponent activity coeffs
      DO IAN = 1, NAN
      DO ICAT = 1, NCAT
         ZBAR           = ( ZP( ICAT ) + ZM( IAN ) ) * 0.5d0
         ZBAR2          = ZBAR * ZBAR
         Y( IAN, ICAT ) = ZBAR2 * AN( IAN ) / I
         X( ICAT, IAN ) = ZBAR2 * CAT( ICAT ) / I
      ENDDO
      ENDDO

      DO IAN = 1, NAN
         F1( IAN ) = 0.0d0
         DO ICAT = 1, NCAT
            F1( IAN ) = F1( IAN ) + X( ICAT, IAN ) * LGAMA0( ICAT, IAN )&
                     + ZOT1 * ZP( ICAT ) * ZM( IAN ) * X( ICAT, IAN )
         ENDDO
      ENDDO

      DO ICAT = 1, NCAT
         F2( ICAT ) = 0.0d0
         DO IAN = 1, NAN
            F2( ICAT ) = F2( ICAT ) + Y( IAN, ICAT ) * LGAMA0(ICAT, IAN)&
                      + ZOT1 * ZP( ICAT ) * ZM( IAN ) * Y( IAN, ICAT )
         ENDDO
      ENDDO

      ! now calculate the multicomponent activity coefficients
      DO IAN  = 1, NAN
      DO ICAT = 1, NCAT

         TA  = -ZOT1 * ZP( ICAT ) * ZM( IAN )
         TB  = ZP( ICAT ) * ZM( IAN ) / ( ZP( ICAT ) + ZM( IAN ) )
         TC  = ( F2( ICAT ) / ZP( ICAT ) + F1( IAN ) / ZM( IAN ) )
         TRM = TA + TB * TC
         
         IF ( TRM .GT. 30.0d0 ) THEN
            GAMA( ICAT, IAN ) = 1.0d+30
         ELSE
            GAMA( ICAT, IAN ) = 10.0d0**TRM
         ENDIF

      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE ACTCOF_new


 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!svn version 2908 for testing if there are significant differences
!will be deleted
      subroutine rpmares_2900 ( SO4, HNO3, NO3, NH3, NH4, RH, TEMP,   &
                           ASO4, ANO3, AH2O, ANH4, GNH3, GNO3,   &
                           ERRMARK,debug_flag) 

!-----------------------------------------------------------------------
!C
!C Description:
!C
!C   ARES calculates the chemical composition of a sulfate/nitrate/
!C   ammonium/water aerosol based on equilibrium thermodynamics.
!C
!C   This code considers two regimes depending upon the molar ratio 
!C   of ammonium to sulfate. 
!C
!C   For values of this ratio less than 2,the code solves a cubic for 
!C   hydrogen ion molality, HPLUS,  and if enough ammonium and liquid
!C   water are present calculates the dissolved nitric acid. For molal
!C   ionic strengths greater than 50, nitrate is assumed not to be present. 
!C   
!C   For values of the molar ratio of 2 or greater, all sulfate is assumed
!C   to be ammonium sulfate and a calculation is made for the presence of
!C   ammonium nitrate.
!C
!C   The Pitzer multicomponent approach is used in subroutine ACTCOF to
!C   obtain the activity coefficients. Abandoned -7/30/97 FSB 

!c   The Bromley method of calculating the activity coefficients is s used
!c    in this version

!c   The calculation of liquid water
!C   is done in subroutine water. Details for both calculations are given
!C   in the respective subroutines.
!C
!C   Based upon MARS due to 
!C   P. Saxena, A.B. Hudischewskyj, C. Seigneur, and J.H. Seinfeld, 
!C   Atmos. Environ., vol. 20, Number 7, Pages 1471-1483, 1986.
!C
!C   and SCAPE due to 
!C   Kim, Seinfeld, and Saxeena, Aerosol Ceience and Technology,
!C   Vol 19, number 2, pages 157-181 and pages 182-198, 1993.
!C
!C NOTE: All concentrations supplied to this subroutine are TOTAL
!C       over gas and aerosol phases
!C
!C Parameters:
!C 
!C  SO4   : Total sulfate in MICROGRAMS/M**3 as sulfate (IN)
!C  HNO3  : Nitric Acid in MICROGRAMS/M**3 as nitric acid (IN)
!C  NO3   : Total nitrate in MICROGRAMS/M**3 as nitric acid (IN)
!C  NH3   : Total ammonia in MICROGRAMS/M**3 as ammonia (IN)
!C  NH4   : Ammonium in MICROGRAMS/M**3 as ammonium (IN)
!C  RH    : Fractional relative humidity (IN)
!C  TEMP  : Temperature in Kelvin (IN)
!C  GNO3  : Gas phase nitric acid in MICROGRAMS/M**3 (OUT)
!C  GNH3  : Gas phase ammonia in MICROGRAMS/M**3 (OUT)
!C  ASO4  : Aerosol phase sulfate in MICROGRAMS/M**3 (OUT) 
!C  ANO3  : Aerosol phase nitrate in MICROGRAMS/M**3 (OUT)
!C  ANH4  : Aerosol phase ammonium in MICROGRAMS/M**3 (OUT)
!C  AH2O  : Aerosol phase water in MICROGRAMS/M**3 (OUT)
!C  NITR  : Number of iterations for obtaining activity coefficients  (OUT) 
!C  NR    : Number of real roots to the cubic in the low ammonia case (OUT)
!C 
!C Revision History:
!C      Who       When        Detailed description of changes
!C   ---------   --------  -------------------------------------------
!C   S.Roselle   11/10/87  Received the first version of the MARS code
!C   S.Roselle   12/30/87  Restructured code
!C   S.Roselle   2/12/88   Made correction to compute liquid-phase 
!C                         concentration of H2O2.  
!C   S.Roselle   5/26/88   Made correction as advised by SAI, for 
!C                         computing H+ concentration.
!C   S.Roselle   3/1/89    Modified to operate with EM2
!C   S.Roselle   5/19/89   Changed the maximum ionic strength from 
!C                         100 to 20, for numerical stability.
!C   F.Binkowski 3/3/91    Incorporate new method for ammonia rich case
!C                         using equations for nitrate budget.
!C   F.Binkowski 6/18/91   New ammonia poor case which 
!C                         omits letovicite.
!C   F.Binkowski 7/25/91   Rearranged entire code, restructured
!C                         ammonia poor case.
!C   F.Binkowski 9/9/91    Reconciled all cases of ASO4 to be output
!C                         as SO4--
!C   F.Binkowski 12/6/91   Changed the ammonia defficient case so that 
!C                         there is only neutralized sulfate (ammonium
!C                         sulfate) and sulfuric acid.
!C   F.Binkowski 3/5/92    Set RH bound on AWAS to 37 % to be in agreement 
!C                          with the Cohen et al. (1987)  maximum molality
!C                          of 36.2 in Table III.( J. Phys Chem (91) page
!C                          4569, and Table IV p 4587.)
!C   F.Binkowski 3/9/92    Redid logic for ammonia defficient case to remove
!C                         possibility for denomenator becoming zero; 
!C                         this involved solving for HPLUS first.
!C                         Note that for a relative humidity
!C                          less than 50%, the model assumes that there is no 
!C                          aerosol nitrate.
!C   F.Binkowski 4/17/95   Code renamed  ARES (AeRosol Equilibrium System)  
!C                          Redid logic as follows
!C                         1. Water algorithm now follows Spann & Richardson
!C                         2. Pitzer Multicomponent method used
!C                         3. Multicomponent practical osmotic coefficient 
!C                            use to close iterations.
!C                         4. The model now assumes that for a water
!C                            mass fraction WFRAC less than 50% there is
!C                            no aerosol nitrate.
!C   F.Binkowski 7/20/95   Changed how nitrate is calculated in ammonia poor 
!C                         case, and changed the WFRAC criterion to 40%.
!C                         For ammonium to sulfate ratio less than 1.0 
!C                         all ammonium is aerosol and no nitrate aerosol
!C                         exists.
!C   F.Binkowski 7/21/95   Changed ammonia-ammonium in ammonia poor case to
!C                         allow gas-phase ammonia to exist. 
!C   F.Binkowski 7/26/95   Changed equilibrium constants to values from 
!C                         Kim et al. (1993) 
!C   F.Binkowski 6/27/96   Changed to new water format
!c   F.Binkowski 7/30/97   Changed to Bromley method for multicomponent 
!c                         activity coefficients. The binary activity coefficients 
!c                         are the same as the previous version
!c   F.Binkowski 8/1/97    Chenged minimum sulfate from 0.0 to 1.0e-6 i.e.
!c                         1 picogram per cubic meter
!C   I.Ackermann 2/23/98   modification for solving precision problems
!C        	  on workstations
!C   I.Ackermann 2/25/99   changed logic as follows:
!c                         If iterations fail, initial values of nitrate
!c                          are retained.
!c                         Total mass budgets are changed to account for gas
!c                         phase returned. (incorporated from FZB's models3
!c                         framework)
!C                         eliminated ratio=5 !! for low to zero sulfate
!C   I.Ackermann 3/17/99   modified ratio = 5 treatment see RB3,p.125
!C
!C-----------------------------------------------------------------------


!...........ARGUMENTS and their descriptions

  real, intent(in) ::  SO4   &     ! Total sulfate in micrograms / m**3 
                      ,HNO3  &     ! Total nitric acid in micrograms / m**3
                      ,NO3   &     ! Total nitrate in micrograms / m**3
                      ,NH3   &     ! Total ammonia in micrograms / m**3
                      ,NH4   &     ! Total ammonium in micrograms / m**3
                      ,RH    &     ! Fractional relative humidity 
                      ,TEMP        ! Temperature in Kelvin 

  real, intent(out)::  ASO4  &     ! Aerosol sulfate in micrograms / m**3 
                      ,ANO3  &     ! Aerosol nitrate in micrograms / m**3
                      ,AH2O  &     ! Aerosol liquid water content water in micrograms / m**3
                      ,ANH4  &     ! Aerosol ammonium in micrograms / m**3
                      ,GNO3  &     ! Gas-phase nitric acid in micrograms / m**3
                      ,GNH3        ! Gas-phase ammonia in micrograms / m**3

  logical, intent(in) :: debug_flag

!C...........INCLUDES and their descriptions
!!      INCLUDE SUBST_CONST          ! constants

!...........PARAMETERS and their descriptions:

      real        MWNACL           ! molecular weight for NaCl
      parameter ( MWNACL = 58.44277 )

!ds moved a bunch upstairs.

      real        MWORG            ! molecular weight for Organic Species
      parameter ( MWORG = 16.0 )

      real        MWCL             ! molecular weight for Chloride  
      parameter ( MWCL = 35.453 )

      real        MWAIR            ! molecular weight for AIR
      parameter ( MWAIR = 28.964 )

      real        MWLCT            ! molecular weight for Letovicite
      parameter ( MWLCT = 3.0 * MWNH4 + 2.0 * MWSO4 + 1.0080 )

      real        MWAS             ! molecular weight for Ammonium Sulfate
      parameter ( MWAS = 2.0 * MWNH4 + MWSO4 )

      real        MWABS            ! molecular weight for Ammonium Bisulfate 
      parameter ( MWABS = MWNH4 + MWSO4 + 1.0080 )


!...........SCRATCH LOCAL VARIABLES and their descriptions:
       
      REAL        fRH              ! Index set to percent relative humidity  
      INTEGER     NITR             ! Number of iterations for activity coefficients
      INTEGER     NNN              ! Loop index for iterations 
      INTEGER     NR               ! Number of roots to cubic equation for HPLUS
      INTEGER     ERRMARK

      real          A0             ! Coefficients and roots of 
      real          A1             ! Coefficients and roots of 
      real          A2             ! Coefficients and roots of 
      real        AA               ! Coefficients and discriminant for quadratic equation for ammonium nitrate
      real        BAL              ! internal variables ( high ammonia case)
      real        BB               ! Coefficients and discriminant for quadratic equation for ammonium nitrate
      real        BHAT             ! Variables used for ammonia solubility 
      real        CC               ! Coefficients and discriminant for quadratic equation for ammonium nitrate
      real        CONVT            ! Factor for conversion of units  
      real        DD               ! Coefficients and discriminant for quadratic equation for ammonium nitrate
      real        DISC             ! Coefficients and discriminant for quadratic equation for ammonium nitrate
      real        EROR             ! Relative error used for convergence test  
      real        FNH3             ! "Free ammonia concentration", that which exceeds TWOSO4       
      real        GAMAAB           ! Activity Coefficient for (NH4+, HSO4-)GAMS( 2,3 )
      real        GAMAAN           ! Activity coefficient for (NH4+, NO3-) GAMS( 2,2 )
      real        GAMAHAT          ! Variables used for ammonia solubility 
      real        GAMANA           ! Activity coefficient for (H+ ,NO3-)   GAMS( 1,2 )
      real        GAMAS1           ! Activity coefficient for (2H+, SO4--) GAMS( 1,1 )
      real        GAMAS2           ! Activity coefficient for (H+, HSO4-)  GAMS( 1,3 )
      real        GAMOLD           ! used for convergence of iteration
      real        GASQD            ! internal variables ( high ammonia case)
      real        HPLUS            ! Hydrogen ion (low ammonia case) (moles / kg water)
      real        K1A              ! Equilibrium constant for ammoniua to ammonium
      real        K2SA             ! Equilibrium constant for sulfate-bisulfate (aqueous)
      real        K3               ! Dissociation constant for ammonium nitrate 
      real        KAN              ! Equilibrium constant for ammonium nitrate (aqueous)
      real        KHAT             ! Variables used for ammonia solubility 
      real        KNA              ! Equilibrium constant for nitric acid (aqueous)   
      real        KPH              ! Henry's Law Constant for ammonia       
      real        KW               ! Equilibrium constant for water dissociation             
      real        KW2              ! Internal variable using KAN 
      real        MAN              ! Nitrate (high ammonia case) (moles / kg water)
      real        MAS              ! Sulfate (high ammonia case) (moles / kg water)
      real        MHSO4            ! Bisulfate (low ammonia case) (moles / kg water)
      real        MNA              ! Nitrate (low ammonia case) (moles / kg water)
      real        MNH4             ! Ammonium (moles / kg water)
      real        MOLNU            ! Total number of moles of all ions
      real        MSO4             ! Sulfate (low ammonia case) (moles / kg water)
      real        PHIBAR           ! Practical osmotic coefficient      
      real        PHIOLD           ! Previous value of practical osmotic coefficient used for convergence of iteration
      real        RATIO            ! Molar ratio of ammonium to sulfate
      real        RK2SA            ! Internal variable using K2SA
      real        RKNA             ! Internal variables using KNA
      real        RKNWET           ! Internal variables using KNA
      real        RR1
      real        RR2
      real        STION            ! Ionic strength
      real        T1               ! Internal variables for temperature corrections
      real        T2               ! Internal variables for temperature corrections
      real        T21              ! Internal variables of convenience (low ammonia case)
      real        T221             ! Internal variables of convenience (low ammonia case)
      real        T3               ! Internal variables for temperature corrections
      real        T4               ! Internal variables for temperature corrections
      real        T6               ! Internal variables for temperature corrections
      real        TNH4             ! Total ammonia and ammonium in micromoles / meter ** 3
      real        TNO3             ! Total nitrate in micromoles / meter ** 3
      real        TOLER1           ! Tolerances for convergence test 
      real        TOLER2           ! Tolerances for convergence test 
      real        TSO4             ! Total sulfate in micromoles / meter ** 3
      real        TWOSO4           ! 2.0 * TSO4  (high ammonia case) (moles / kg water)
      real        WFRAC            ! Water mass fraction 
      real        WH2O             ! Aerosol liquid water content (internally)
                                   ! micrograms / meter **3 on output
                                   ! internally it is 10 ** (-6) kg (water) / meter ** 3
                                   ! the conversion factor (1000 g = 1 kg) is applied 
                                   ! for AH2O output
      real        WSQD             ! internal variables ( high ammonia case)
      real        XNO3             ! Nitrate aerosol concentration in micromoles / meter ** 3
      real        XXQ              ! Variable used in quadratic solution
      real        YNH4             ! Ammonium aerosol concentration in micromoles / meter** 3
      real        ZH2O             ! Water variable saved in case ionic strength too high.
      real        ZSO4             ! Total sulfate molality - mso4 + mhso4 (low ammonia case) (moles / kg water)

      real        CAT( 2 )         ! Array for cations (1, H+); (2, NH4+) (moles / kg water)
      real        AN ( 3 )         ! Array for anions (1, SO4--); (2, NO3-); (3, HSO4-)  (moles / kg water) 
      real        CRUTES( 3 )      ! Coefficients and roots of 
      real        GAMS( 2, 3 )     ! Array of activity coefficients 
      real        MINSO4           ! Minimum value of sulfate laerosol concentration
       parameter( MINSO4 = 1.0E-6 / MWSO4 ) 
      real        MINNO3
       parameter( MINNO3 = 1.0E-6 / MWNO3 )    !2/25/99 IJA
!st      real        FLOOR
!st       parameter( FLOOR = 1.0E-30) ! minimum concentration       
!2/25/99 IJA
! FSB New variables Total ammonia and nitrate mass concentrations
      real  TMASSNH3  ! Total ammonia (gas and particle)
                      !  as ammonia gas [ug /m**3]
      real  TMASSHNO3 ! Total nitrate (vapor and particle) as
                      !  nitric acid vapor [ug /m**3]

!for RATIO between RATIO_Low and RATIO_High, consider the species in two boxes one with
! RATIO_Low and the other with  RATIO_High. Then take a linear combination of the results in each box
      logical ::DO_RATIO_High_2,DO_RATIO_Low_2
      real, parameter ::RATIO_Low = 1.95,RATIO_High= 2.05 !somewhat arbitrarily
      real ::ASO4_High, ANO3_High, AH2O_High, ANH4_High, GNH3_High, GNO3_High
      real ::ASO4_High1, ANO3_High1, AH2O_High1, ANH4_High1, GNH3_High1, GNO3_High1,diff1
      real ::ASO4_Low, ANO3_Low, AH2O_Low, ANH4_Low, GNH3_Low, GNO3_Low
      real ::ASO4_Low1, ANO3_Low1, AH2O_Low1, ANH4_Low1, GNH3_Low1, GNO3_Low1,diff2
      real::TSO4_HighA,TSO4_LowA,High_Factor,X
!-----------------------------------------------------------------------
!  begin body of subroutine RPMARES
                                                                         
!ASO4=FLOOR;ANO3=FLOOR;AH2O=FLOOR;ANH4=FLOOR;GNO3=FLOOR;GNH3=FLOOR 
!Initialise the output variables

      ASO4=0.0;ANO3=0.0;AH2O=0.0;ANH4=0.0;GNO3=0.0;GNH3=0.0 

      ASO4 = SO4   !ds from RPM

!...convert into micromoles/m**3
 
!..iamodels3 merge NH3/NH4 , HNO3,NO3 here
      TSO4 = MAX( 0.0, SO4 / MWSO4  )
      TNO3 = MAX( 0.0, (NO3 / MWNO3 + HNO3 / MWHNO3) )
      TNH4 = MAX( 0.0, (NH3 / MWNH3 + NH4 / MWNH4)  )

!2/25/99 IJA
!      TMASSNH3  = MAX(0.0, NH3 + (MWNH3 / MWNH4)  * NH4 )
!      TMASSHNO3 = MAX(0.0, NO3 + (MWHNO3 / MWNO3) * NO3 )

      TMASSNH3  = MAX(0.0, NH3 +  NH4 )
      TMASSHNO3 = MAX(0.0, HNO3 + NO3 )
 
!...now set humidity index fRH as a percent

!st      IRH = NINT( 100.0 * RH )
         fRH = RH 
!...Check for valid fRH

       fRH = MAX( 0.01, fRH )
       fRH = MIN( 0.99, fRH )

!...Specify the equilibrium constants at  correct
!...  temperature.  Also change units from ATM to MICROMOLE/M**3 (for KAN,
!...  KPH, and K3 )
!...  Values from Kim et al. (1993) except as noted.
 
      CONVT = 1.0 / ( 0.082 * TEMP ) 
      T6 = 0.082E-9 *  TEMP
      T1 = 298.0 / TEMP
      T2 = ALOG( T1 )
      T3 = T1 - 1.0
      T4 = 1.0 + T2 - T1
      KNA  = 2.511E+06 *  EXP(  29.17 * T3 + 16.83 * T4 ) * T6
      K1A  = 1.805E-05 *  EXP(  -1.50 * T3 + 26.92 * T4 )
      K2SA = 1.015E-02 *  EXP(   8.85 * T3 + 25.14 * T4 )
      KW   = 1.010E-14 *  EXP( -22.52 * T3 + 26.92 * T4 )
      KPH  = 57.639    *  EXP(  13.79 * T3 - 5.39  * T4 ) * T6
!!!      K3   =  5.746E-17 * EXP( -74.38 * T3 + 6.12  * T4 ) * T6 * T6  
      KHAT =  KPH * K1A / KW  
      KAN  =  KNA * KHAT  

!...Compute temperature dependent equilibrium constant for NH4NO3
!...  ( from Mozurkewich, 1993)

      K3 = EXP( 118.87  - 24084.0 / TEMP -  6.025  * ALOG( TEMP ) )

!...Convert to (micromoles/m**3) **2

      K3 = K3 * CONVT * CONVT

      WH2O   = 0.0
      STION  = 0.0
      AH2O   = 0.0
      MAS    = 0.0
      MAN    = 0.0
      HPLUS  = 0.0
      TOLER1 = 0.00001
      TOLER2 = 0.001
      NITR   = 0
      NR     = 0
      RATIO  = 0.0
      GAMAAN = 1.0
      GAMOLD = 1.0

!ds from RPM, but removed FLOOR : (slighty different logic)
      if( (TSO4 < MINSO4 ) .and. (TNO3 < MINNO3) ) then
        !print *, "DSX  HERE"
          ASO4 = SO4   ! MAX..
          ANO3 = NO3   ! MAX..      
          WH2O = 0.0
          AH2O = 0.0
          GNH3 = NH3   ! MAX(FLOOR,NH3)
          GNO3 = NO3   ! MAX(FLOOR,NO3)
          RETURN
       END IF
!ds end rpm


!...set the ratio according to the amount of sulfate and nitrate

      IF ( TSO4 > MINSO4 ) THEN
        RATIO = TNH4 / TSO4

!...If there is no sulfate and no nitrate, there can be no ammonium
!...  under the current paradigm. Organics are ignored in this version.

      ELSE 
      
        IF ( TNO3 <= MINNO3 ) THEN

! *** If there is very little sulfate and nitrate set concentrations
!      to a very small value and return.    
! Jun 2012, Note these values are set in the initialisation
          ASO4 = SO4 ! MAX(FLOOR, ASO4)
          ANO3 = NO3 ! MAX(FLOOR, ANO3 )          
          WH2O = 0.0
          AH2O = 0.0
          GNH3 = NH3 ! MAX(FLOOR,GNH3)
          GNO3 = NO3 ! MAX(FLOOR,GNO3)
          RETURN
       END IF
       
!...For the case of no sulfate and nonzero nitrate, set ratio to 5
!...  to send the code to the high ammonia case if there is more
!...  ammonia than sulfate, otherwise send to low ammonia case.

       IF (TNH4 > TSO4) THEN
         RATIO = 5.0        !this is a high ammonia case with low sulfur
       ELSE
         RATIO = 1.        !this is a low ammonia case with low sulfur
       END IF
 
      END IF 


      DO_RATIO_High_2=.false.
      DO_RATIO_Low_2=.false.
      ASO4_High = 0.0
      ANO3_High = 0.0
      AH2O_High = 0.0
      ANH4_High = 0.0
      GNH3_High = 0.0
      GNO3_High = 0.0
      ASO4_Low = 0.0
      ANO3_Low = 0.0
      AH2O_Low = 0.0
      ANH4_Low = 0.0
      GNH3_Low = 0.0
      GNO3_Low = 0.0
      TSO4_HighA = TSO4
      TSO4_LowA  = TSO4
       IF ( RATIO > 2.0 )then
         DO_RATIO_High_2=.true.
         High_Factor=1.0
      else
         DO_RATIO_Low_2=.true.
         High_Factor=0.0
      end if

      IF ( RATIO >RATIO_Low  .and.  RATIO < RATIO_High .and. MARS_RATIO_SMOOTH)then
         DO_RATIO_High_2=.true.
         DO_RATIO_Low_2=.true.
         TSO4_HighA=TSO4*Ratio/RATIO_High
         TSO4_LowA=TSO4*Ratio/RATIO_Low
!         High_Factor=(RATIO-RATIO_Low)/(RATIO_High-RATIO_Low)
         High_Factor=(RATIO*RATIO_High-RATIO_Low*RATIO_High)/(RATIO*RATIO_High-RATIO*RATIO_Low)
      end if

!....................................
!......... High Ammonia Case ........
!....................................
 
!      IF ( RATIO > 2.0 ) THEN    ! NH4/SO4 > 2
      IF ( DO_RATIO_High_2 ) THEN    ! NH4/SO4 > 2
        GAMAAN = 0.1

!...Set up twice the sulfate for future use.

        TWOSO4 = 2.0 * TSO4_HighA
        XNO3 = 0.0            
        YNH4 = TWOSO4

!...Treat different regimes of relative humidity 

!...ZSR relationship is used to set water levels. Units are
!...  10**(-6) kg water/ (cubic meter of air)
!...  start with ammomium sulfate solution without nitrate

      CALL awater(fRH,TSO4_HighA,YNH4,TNO3,AH2O ) !**** note TNO3
        WH2O = 1.0E-3 * AH2O  
        ASO4_High = TSO4_HighA * MWSO4
        ANO3_High = 0.0
        ANH4_High = YNH4 * MWNH4
if(debug_flag) call datewrite("MARS debug ", -1,(/ ASO4_High, ANH4_High, AH2O /) )

if ( DEBUG%EQUIB ) then
  if( ASO4_High + ANH4_High +  AH2O < 1.0-10 ) then
     call datewrite("MARS failing? ", -1,(/ ASO4_High, ANH4_High, AH2O /) )
     print *, "MARS PROB ", ASO4_High, ANH4_High, AH2O, TSO4_HighA, YNH4
     call CheckStop("MARS")
  end if
end if
        WFRAC = (AH2O + FLOOR)  / ( ASO4_High + ANH4_High +  AH2O + FLOOR  )
        !ds WFRAC = AH2O / ( ASO4 + ANH4 +  AH2O + FLOOR  )
        !CRUDE FIX? WFRAC = AH2O / ( ASO4 + ANH4 +  AH2O )
!!!!       IF ( WFRAC == 0.0 )  RETURN   ! No water       
        IF ( WFRAC < 0.2 ) THEN
 
!..."dry" ammonium sulfate and ammonium nitrate
!...  compute free ammonia

          FNH3 = TNH4 - TWOSO4
          CC = TNO3 * FNH3 - K3

!...check for not enough to support aerosol      

          !dsjIF ( CC <= 0.0 ) THEN
          IF ( CC < FLOOR ) THEN
            XNO3 = 0.0
          ELSE
            AA = 1.0
            BB = -( TNO3 + FNH3 ) 
            DISC = BB * BB - 4.0 * CC

!...Check for complex roots of the quadratic
!...  set nitrate to zero and RETURN if complex roots are found
!2/25/99 IJA

            !DS IF ( DISC < 0.0 ) THEN
            !dsDSdIF ( DISC < FLOOR ) THEN
            IF ( DISC < 0.0 ) THEN
if( DEBUG%EQUIB .and. debug_flag ) print *, "MARS DISC NEG ", XNO3, WH2O, DISC
              XNO3 = 0.0
              AH2O_High = 1000.0 * WH2O
              YNH4 = TWOSO4
              GNO3_High = HNO3
              ASO4_High = TSO4_HighA * MWSO4
              ANO3_High = NO3
              ANH4_High = YNH4 * MWNH4
              GNH3_High = TMASSNH3 - ANH4
              if( GNH3 < 0.0 ) then
                 print *, " NEG GNH3", TWOSO4, ANH4, TMASSNH3
                 call CheckStop("NEG GNH3")
              end if
!              RETURN
          goto 333
            END IF

!...to get here, BB .lt. 0.0, CC .gt. 0.0 always      

            DD = SQRT( DISC )
            XXQ = -0.5 * ( BB + SIGN ( 1.0, BB ) * DD )

!...Since both roots are positive, select smaller root.      

            XNO3 = MIN( XXQ / AA, CC / XXQ )
          
          END IF
!2/25/99 IJA
          AH2O_High = 1000.0 * WH2O
          YNH4 = TWOSO4 + XNO3
          ASO4_High = TSO4_HighA * MWSO4
          !dsSAFE ANO3 = XNO3 * MWNO3
          ANO3_High = min(XNO3 * MWNO3, TMASSHNO3 )
          !dsSAFE ANH4 = YNH4 * MWNH4 ! ds should be safe as NH4/SO4 >2, but anyway:
          ANH4_High = min(YNH4 * MWNH4, TMASSNH3 )  ! ds should be safe as NH4/SO4 >2, but anyway:
          GNH3_High = TMASSNH3 - ANH4_High
          GNO3_High = TMASSHNO3 - ANO3_High
          !    if( GNH3 < 0.0 .or. GNO3 < 0.0 ) then
          !       print *, " NEG GNH3 GNO3", TWOSO4, ANH4, TMASSNH3, ANO3, TMASSHNO3
          !       call CheckStop("NEG GNH3 GNO3")
          !    end if
!          RETURN
          goto 333
        END IF

!...liquid phase containing completely neutralized sulfate and
!...  some nitrate.  Solve for composition and quantity.
 
        MAS = TSO4_HighA / max(WH2O,FLOOR)
        MAN = 0.0
        XNO3 = 0.0
        YNH4 = TWOSO4
        PHIOLD = 1.0

!...Start loop for iteration
 
!...The assumption here is that all sulfate is ammonium sulfate,
!...  and is supersaturated at lower relative humidities.
        diff1=1.0/TOLER1

        DO 1501 NNN = 1, 150 
          NITR = NNN
          GASQD = GAMAAN * GAMAAN
          WSQD = WH2O * WH2O
          KW2 = KAN * WSQD / GASQD
          AA = 1.0 - KW2
          BB = TWOSO4 + KW2 * ( TNO3 + TNH4 - TWOSO4 )
          CC = -KW2 * TNO3 * ( TNH4 - TWOSO4 )

!...This is a quadratic for XNO3 [MICROMOLES / M**3] of nitrate in solution.

          DISC = BB * BB - 4.0 * AA * CC

          if( DEBUG%EQUIB ) then
            MAXNNN1 = NNN
            !if( MAXNNN1 > 140) print "(a,i4,9es12.3)", "NNN1 ", NNN, DISC, TNO3, TNH4, TWOSO4
          end if
!...Check for complex roots, retain inital values and RETURN
!2/25/99 IJA

          !DS IF ( DISC < 0.0 ) THEN
          !dsDS IF ( DISC < FLOOR ) THEN
          IF ( DISC < 0.0 ) THEN
if( DEBUG%EQUIB .and. debug_flag ) print *, "MARS DISC NEG2 ", XNO3, WH2O, DISC
            XNO3 = 0.0
            AH2O_High = 1000.0 * WH2O
            YNH4 = TWOSO4
            GNO3_High = HNO3
            ASO4_High = TSO4_HighA * MWSO4
            ANO3_High = NO3
            !ANH4 = YNH4 * MWNH4
            ANH4_High = min( YNH4 * MWNH4, TMASSNH3)  ! ds added "min"
            GNH3_High = TMASSNH3 - ANH4_High

              !WRITE( 10, * ) ' COMPLEX ROOTS '
!            RETURN
                goto 333
          END IF
! 2/25/99 IJA

! Deal with degenerate case (yoj)

          !DS IF ( AA /= 0.0 ) THEN
          IF ( abs(AA) > FLOOR  ) THEN
             if( DEBUG%EQUIB .and. debug_flag ) print "(a,9es11.3)", "MARS DEGEN  ",  XNO3, WH2O, DISC, AA, BB, CC
             DD = SQRT( DISC )
             XXQ = -0.5 * ( BB + SIGN ( 1.0, BB ) * DD )

             RR1 = XXQ / AA
             RR2 = CC / XXQ

             
             !...choose minimum positve root         
             
             IF ( ( RR1 * RR2 ) < 0.0 ) THEN
                if( DEBUG%EQUIB .and. debug_flag ) print "(a,10es10.3)", "MARS RR1*RR2  ", XNO3, WH2O, DISC, RR1, RR2
                XNO3 = MAX( RR1, RR2 )
             ELSE if(MIN( RR1, RR2 )>0.0)then
                XNO3 = MIN( RR1, RR2 )
             ELSE!two negative roots !DS PW added 4th July 2012

                           !--------------------- return copied from above
                XNO3 = 0.0
                AH2O_High = 1000.0 * WH2O
                YNH4 = TWOSO4
                GNO3_High = HNO3
                ASO4_High = TSO4_HighA * MWSO4
                ANO3_High = NO3
                !ds ANH4 = YNH4 * MWNH4
                ANH4_High = min( YNH4 * MWNH4, TMASSNH3)  ! ds added "min"
                GNH3_High = TMASSNH3 - ANH4_High
                if( DEBUG%EQUIB .and. debug_flag ) WRITE( *, * ) ' TWO NEG ROOTS '
!                RETURN
                goto 333
               
             END IF
          ELSE
             XNO3 = - CC / BB   ! AA equals zero here
if( DEBUG%EQUIB .and. debug_flag ) print "(a,4es10.3)", "MARS NONDEGEN  ",  AA, BB, CC, XNO3
          END IF

          XNO3 = MIN( XNO3, TNO3 )

!...This version assumes no solid sulfate forms (supersaturated ) 
!...  Now update water

          CALL AWATER ( fRH, TSO4_HighA, YNH4, XNO3, AH2O)

!...ZSR relationship is used to set water levels. Units are
!...  10**(-6) kg water/ (cubic meter of air)
!...  The conversion from micromoles to moles is done by the units of WH2O.

          WH2O = 1.0E-3 * AH2O

!...Ionic balance determines the ammonium in solution.

          MAN = XNO3 / max(WH2O,FLOOR)
          MAS = TSO4_HighA / max(WH2O,FLOOR)
          MNH4 = 2.0 * MAS + MAN
          YNH4 = MNH4 * WH2O

 !st ...  FIXING
   if(MNH4<0. .or. MAS<0. .or. MAN<0.) then
      MNH4 = 1.e-30
      MAS  = 1.e-30
      MAN  = 1.e-30
   end if

!...MAS, MAN and MNH4 are the aqueous concentrations of sulfate, nitrate,
!...  and ammonium in molal units (moles/(kg water) ).

          STION = 3.0 * MAS + MAN
          CAT( 1 ) = 0.0
          CAT( 2 ) = MNH4 
          AN ( 1 ) = MAS
          AN ( 2 ) = MAN
          AN ( 3 ) = 0.0
          CALL ACTCOF ( CAT, AN, GAMS, MOLNU, PHIBAR , ERRMARK,1,debug_flag)
          GAMAAN = GAMS( 2, 2 )

!...Use GAMAAN for convergence control

          EROR = ABS( GAMOLD - GAMAAN ) / GAMOLD
          GAMOLD = GAMAAN

!...Check to see if we have a solution

          IF ( EROR <= TOLER1 ) THEN 
!!!            WRITE( 11, * ) RH, STION, GAMS( 1, 1 ),GAMS( 1, 2 ), GAMS( 1, 3 ),
!!!     &      GAMS( 2, 1 ), GAMS( 2, 2 ), GAMS( 2, 3 ), PHIBAR
! 2/25/99 IJA

            ASO4_High = TSO4_HighA * MWSO4
            ANO3_High = min(XNO3 * MWNO3,TMASSHNO3)
            ANH4_High = min( YNH4 * MWNH4, TMASSNH3 )
            GNO3_High = TMASSHNO3  - ANO3_High
            GNH3_High = TMASSNH3   - ANH4_High
            AH2O_High = 1000.0 * WH2O

!            RETURN
            goto 333
          END IF


1501    CONTINUE

!...If after NITR iterations no solution is found, then:
! FSB retain the initial values of nitrate particle and vapor
! 2/25/99 IJA
        ASO4_High = TSO4_HighA * MWSO4
        ANO3_High = NO3
        XNO3 = NO3 / MWNO3
        YNH4 = TWOSO4
!        ANH4 = YNH4 * MWNH4
        ANH4_High = min( YNH4 * MWNH4, TMASSNH3 )  ! ds pw added "min"
        CALL AWATER ( fRH, TSO4_HighA, YNH4, XNO3, AH2O_High)
        GNO3_High = HNO3
        GNH3_High = TMASSNH3 - ANH4_High
!        RETURN
        goto 333

333 CONTINUE !finished DO_RATIO_High_2

!      ELSE
     ENDIF
      IF ( DO_RATIO_Low_2 ) THEN    ! NH4/SO4 <= 2
       
!......................................
!......... Low Ammonia Case ...........
!......................................
      
!...coded by Dr. Francis S. Binkowski 12/8/91.(4/26/95)
      ! modified 8/28/98

!...All cases covered by this logic
 
        WH2O = 0.0
        CALL AWATER ( fRH, TSO4_LowA, TNH4, TNO3, AH2O )
        WH2O = 1.0E-3 * AH2O
        ZH2O = AH2O

!...convert 10**(-6) kg water/(cubic meter of air) to micrograms of water
!...  per cubic meter of air (1000 g = 1 kg)
! 2/25/99 IJA 
        ASO4_Low = TSO4_LowA * MWSO4
        ANH4_Low = TNH4 * MWNH4
        !dsSAFE ANO3 = NO3
        ANO3_Low = min( NO3, TMASSHNO3 )
        GNO3_Low = TMASSHNO3 - ANO3_Low
        GNH3_Low = FLOOR
        AH2O_Low = 1.0E3 *WH2O

!...Check for zero water.      

        !ds IF ( WH2O == 0.0 ) RETURN
        IF ( abs(WH2O) < FLOOR ) goto 111!RETURN
        ZSO4 = TSO4_LowA / max(WH2O,FLOOR)

!...ZSO4 is the molality of total sulfate i.e. MSO4 + MHSO4      

!!!         IF ( ZSO4 > 11.0 ) THEN

!...do not solve for aerosol nitrate for total sulfate molality
!...  greater than 11.0 because the model parameters break down
!...  greater than  9.0 because the model parameters break down

        IF ( ZSO4 > 9.0 ) THEN   ! 18 June 97
           goto 111
          !RETURN
        END IF

!...First solve with activity coeffs of 1.0, then iterate.

        PHIOLD = 1.0
        GAMANA = 1.0
        GAMAS1 = 1.0
        GAMAS2 = 1.0
        GAMAAB = 1.0
        GAMOLD = 1.0

!...All ammonia is considered to be aerosol ammonium. 

        MNH4 = TNH4 / max(WH2O,FLOOR)

!...MNH4 is the molality of ammonium ion.

        YNH4 = TNH4
      
!...loop for iteration
 
        DO 1601 NNN = 1, 150
          if( DEBUG%EQUIB ) then
            if(NNN > MAXNNN2 ) MAXNNN2 = NNN
            !if( MAXNNN2 > 140) print *, "NNN2 ", NNN, TNO3, TNH4, TWOSO4
          end if
          NITR = NNN

!...set up equilibrium constants including activities
!...  solve the system for hplus first then sulfate & nitrate

          RK2SA = K2SA * GAMAS2 * GAMAS2 / ( GAMAS1 * GAMAS1 * GAMAS1 )
          RKNA = KNA / ( GAMANA * GAMANA )
          RKNWET = RKNA * WH2O       
          T21  = ZSO4 - MNH4
          T221 = ZSO4 + T21

!...set up coefficients for cubic       

          A2 = RK2SA + RKNWET - T21
          A1 = RK2SA * RKNWET - T21 * ( RK2SA + RKNWET )    &
             - RK2SA * ZSO4 - RKNA * TNO3
          A0 = - (T21 * RK2SA * RKNWET                      &
             + RK2SA * RKNWET * ZSO4 + RK2SA * RKNA * TNO3 )   
         
          CALL CUBIC ( A2, A1, A0, NR, CRUTES )
       
!...Code assumes the smallest positive root is in CRUTES(1)
 
          HPLUS = CRUTES( 1 )
          BAL = HPLUS **3 + A2 * HPLUS**2 + A1 * HPLUS + A0
          MSO4 = RK2SA * ZSO4 / ( HPLUS + RK2SA )   ! molality of sulfate ion
          MHSO4 = ZSO4 - MSO4                       ! molality of bisulfate ion
          MNA = RKNA * TNO3 / ( HPLUS + RKNWET )    ! molality of nitrate ion
          MNA = MAX( 0.0, MNA )
          MNA = MIN( MNA, TNO3 / max(WH2O,FLOOR) )
          XNO3 = MNA * WH2O
          !ds ANO3 = MNA * WH2O * MWNO3
          ANO3_Low = min( TMASSHNO3, MNA * WH2O * MWNO3)
! 2/25/99 IJA
          GNO3_Low = TMASSHNO3 - ANO3_Low
          if( DEBUG%EQUIB ) then
             if (GNO3_Low < 0.0 ) call CheckStop("NNN2 GNO3 NEG")
          end if
        
!...Calculate ionic strength      

          STION = 0.5 * ( HPLUS + MNA + MNH4 + MHSO4 + 4.0 * MSO4 )
          
!...Update water

          CALL AWATER ( fRH, TSO4_LowA, YNH4, XNO3, AH2O )

!...Convert 10**(-6) kg water/(cubic meter of air) to micrograms of water
!...  per cubic meter of air (1000 g = 1 kg)                       

          WH2O = 1.0E-3 * AH2O 
          CAT( 1 ) = HPLUS
          CAT( 2 ) = MNH4
          AN ( 1 ) = MSO4
          AN ( 2 ) = MNA
          AN ( 3 ) = MHSO4

          CALL ACTCOF ( CAT, AN, GAMS, MOLNU, PHIBAR, ERRMARK,2,debug_flag)

          GAMANA = GAMS( 1, 2 )
          GAMAS1 = GAMS( 1, 1 )
          GAMAS2 = GAMS( 1, 3 )
          GAMAAN = GAMS( 2, 2 )

          GAMAHAT = ( GAMAS2 * GAMAS2 / ( GAMAAB * GAMAAB ) )
          BHAT = KHAT * GAMAHAT 
!!!          EROR = ABS ( ( PHIOLD - PHIBAR ) / PHIOLD )
!!!          PHIOLD = PHIBAR
          EROR = ABS ( GAMOLD - GAMAHAT ) / GAMOLD 
          GAMOLD = GAMAHAT   

!...write out molalities and activity coefficient
!...  and return with good solution
          AH2O_Low = 1.0E3 *WH2O

          IF ( EROR <= TOLER2 ) THEN
!!!            WRITE(12,*) RH, STION,HPLUS,ZSO4,MSO4,MHSO4,MNH4,MNA
!!!            WRITE(11,*) RH, STION, GAMS(1,1),GAMS(1,2),GAMS(1,3),
!!!     &                  GAMS(2,1),GAMS(2,2),GAMS(2,3), PHIBAR

           goto 111
            !RETURN
          END IF

1601    CONTINUE     

!...after NITR iterations, failure to solve the system, no ANO3
! 2/25/99 IJA
        ANH4_Low = TNH4 * MWNH4
        GNH3_Low = FLOOR
        GNO3_Low = HNO3
        ANO3_Low = NO3
        CALL AWATER ( fRH, TSO4_LowA, TNH4, TNO3, AH2O_Low )      
!        RETURN
           goto 111
111        continue
   
      END IF   ! ratio .gt. 2.0

      ASO4 = ASO4_High*High_Factor + ASO4_Low*(1.0-High_Factor)
      ANO3 = ANO3_High*High_Factor + ANO3_Low*(1.0-High_Factor)
      AH2O = AH2O_High*High_Factor + AH2O_Low*(1.0-High_Factor)
      ANH4 = ANH4_High*High_Factor + ANH4_Low*(1.0-High_Factor)
      GNH3 = GNH3_High*High_Factor + GNH3_Low*(1.0-High_Factor)
      GNO3 = GNO3_High*High_Factor + GNO3_Low*(1.0-High_Factor)

      end subroutine rpmares_2900 ! end RPMares

!>-------------------------------------------------------------------------------<


 end module MARS_mod

