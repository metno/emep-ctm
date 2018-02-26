
!-> SO4 

      P =   &
          2e-12*AQRCK(ICLOHSO2,K) * xnew(OH ) * xnew(SO2 )  &
        + AQRCK(ICLRC1,K) * xnew(SO2 ) * xnew(H2O2 )  &
        + AQRCK(ICLRC2,K) * xnew(SO2 ) * xnew(O3 )  &
        + AQRCK(ICLRC3,K) * xnew(SO2 )  &
        + rcemis(SO4,k) 
      ! L = 0.0
 

      xnew(SO4)=  xold(SO4) + dt2 * P

!-> NH3 

      P =   &
          rcemis(NH3,k)  &
        + 0 !Skip bio rate since rcemis exists 
      ! L = 0.0
 

      xnew(NH3)=  xold(NH3) + dt2 * P

!-> NO3_F 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for NO3_F! xnew(NO3_F)=  max(0.0, xold(NO3_F))

!-> NO3_C 

      P =   &
          rct(80,k) * xnew(HNO3 )  &
        + rct(81,k) * xnew(HNO3 ) 
      ! L = 0.0
 

      xnew(NO3_C)=  xold(NO3_C) + dt2 * P

!-> NH4_F 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for NH4_F! xnew(NH4_F)=  max(0.0, xold(NH4_F))

!-> DUMMY 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for DUMMY! xnew(DUMMY)=  max(0.0, xold(DUMMY))

!-> ASH_F 

      P =   &
          rcemis(ASH_F,k) 
      ! L = 0.0
 

      xnew(ASH_F)=  xold(ASH_F) + dt2 * P

!-> ASH_C 

      P =   &
          rcemis(ASH_C,k) 
      ! L = 0.0
 

      xnew(ASH_C)=  xold(ASH_C) + dt2 * P

!-> FFIRE_CO 

      P =   &
          rcemis(FFIRE_CO,k) 

      L =  &
          rct(24,k)*xnew(OH) 

      xnew(FFIRE_CO)=  ( xold(FFIRE_CO) + dt2 * P) /(1.0 + dt2*L )

!-> POM_F_WOOD 

      P =   &
          rcemis(POM_F_WOOD,k) 
      ! L = 0.0
 

      xnew(POM_F_WOOD)=  xold(POM_F_WOOD) + dt2 * P

!-> POM_F_FFUEL 

      P =   &
          rcemis(POM_F_FFUEL,k) 
      ! L = 0.0
 

      xnew(POM_F_FFUEL)=  xold(POM_F_FFUEL) + dt2 * P

!-> POM_C_FFUEL 

      P =   &
          rcemis(POM_C_FFUEL,k) 
      ! L = 0.0
 

      xnew(POM_C_FFUEL)=  xold(POM_C_FFUEL) + dt2 * P

!-> EC_F_WOOD_NEW 

      P =   &
          rcemis(EC_F_WOOD_NEW,k) 

      L =  &
          rct(110,k) 

      xnew(EC_F_WOOD_NEW)=  ( xold(EC_F_WOOD_NEW) + dt2 * P) /(1.0 + dt2*L )

!-> EC_F_WOOD_AGE 

      P =   &
          rcemis(EC_F_WOOD_AGE,k)  &
        + rct(110,k) * xnew(EC_F_WOOD_NEW ) 
      ! L = 0.0
 

      xnew(EC_F_WOOD_AGE)=  xold(EC_F_WOOD_AGE) + dt2 * P

!-> EC_C_WOOD 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for EC_C_WOOD! xnew(EC_C_WOOD)=  max(0.0, xold(EC_C_WOOD))

!-> EC_F_FFUEL_NEW 

      P =   &
          rcemis(EC_F_FFUEL_NEW,k) 

      L =  &
          rct(111,k) 

      xnew(EC_F_FFUEL_NEW)=  ( xold(EC_F_FFUEL_NEW) + dt2 * P) /(1.0 + dt2*L )

!-> EC_F_FFUEL_AGE 

      P =   &
          rcemis(EC_F_FFUEL_AGE,k)  &
        + rct(111,k) * xnew(EC_F_FFUEL_NEW ) 
      ! L = 0.0
 

      xnew(EC_F_FFUEL_AGE)=  xold(EC_F_FFUEL_AGE) + dt2 * P

!-> EC_C_FFUEL 

      P =   &
          rcemis(EC_C_FFUEL,k) 
      ! L = 0.0
 

      xnew(EC_C_FFUEL)=  xold(EC_C_FFUEL) + dt2 * P

!-> REMPPM25 

      P =   &
          rcemis(REMPPM25,k) 
      ! L = 0.0
 

      xnew(REMPPM25)=  xold(REMPPM25) + dt2 * P

!-> REMPPM_C 

      P =   &
          rcemis(REMPPM_C,k) 
      ! L = 0.0
 

      xnew(REMPPM_C)=  xold(REMPPM_C) + dt2 * P

!-> FFIRE_OM 

      P =   &
          rcemis(FFIRE_OM,k) 
      ! L = 0.0
 

      xnew(FFIRE_OM)=  xold(FFIRE_OM) + dt2 * P

!-> FFIRE_BC 

      P =   &
          rcemis(FFIRE_BC,k) 
      ! L = 0.0
 

      xnew(FFIRE_BC)=  xold(FFIRE_BC) + dt2 * P

!-> FFIRE_REMPPM25 

      P =   &
          rcemis(FFIRE_REMPPM25,k) 
      ! L = 0.0
 

      xnew(FFIRE_REMPPM25)=  xold(FFIRE_REMPPM25) + dt2 * P

!-> OM25_BGND 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for OM25_BGND! xnew(OM25_BGND)=  max(0.0, xold(OM25_BGND))

!-> OM25_P 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for OM25_P! xnew(OM25_P)=  max(0.0, xold(OM25_P))

!-> SQT_SOA_NV 

      P =   &
          0.01534*rcemis(BIOTERP,k) 
      ! L = 0.0
 

      xnew(SQT_SOA_NV)=  xold(SQT_SOA_NV) + dt2 * P

!-> ASOC_NG100 

      P =   &
          rct(94,k)*xnew(OH) * xnew(ASOC_UG1 ) 
      ! L = 0.0
 

      xnew(ASOC_NG100)=  xold(ASOC_NG100) + dt2 * P

!-> ASOC_UG1 

      P =   &
          YCOXY(0) * 1.36e-11*xnew(OXYL)*xnew(OH)  &
        + YCALK(0) * rct(43,k)*xnew(OH)*xnew(NC4H10)  &
        + YCOLE(0) * rct(77,k)*xnew(OH)*xnew(C3H6)  &
        + rct(95,k)*xnew(OH) * xnew(ASOC_UG10 ) 

      L =  &
          rct(94,k)*xnew(OH) 

      xnew(ASOC_UG1)=  ( xold(ASOC_UG1) + dt2 * P) /(1.0 + dt2*L )

!-> ASOC_UG10 

      P =   &
          YCOXY(1) * 1.36e-11*xnew(OXYL)*xnew(OH)  &
        + YCALK(1) * rct(43,k)*xnew(OH)*xnew(NC4H10)  &
        + YCOLE(1) * rct(77,k)*xnew(OH)*xnew(C3H6)  &
        + rct(96,k)*xnew(OH) * xnew(ASOC_UG1E2 ) 

      L =  &
          rct(95,k)*xnew(OH) 

      xnew(ASOC_UG10)=  ( xold(ASOC_UG10) + dt2 * P) /(1.0 + dt2*L )

!-> ASOC_UG1E2 

      P =   &
          YCOXY(2) * 1.36e-11*xnew(OXYL)*xnew(OH)  &
        + YCALK(2) * rct(43,k)*xnew(OH)*xnew(NC4H10)  &
        + YCOLE(2) * rct(77,k)*xnew(OH)*xnew(C3H6)  &
        + rct(97,k)*xnew(OH) * xnew(ASOC_UG1E3 ) 

      L =  &
          rct(96,k)*xnew(OH) 

      xnew(ASOC_UG1E2)=  ( xold(ASOC_UG1E2) + dt2 * P) /(1.0 + dt2*L )

!-> ASOC_UG1E3 

      P =   &
          YCOXY(3) * 1.36e-11*xnew(OXYL)*xnew(OH)  &
        + YCALK(3) * rct(43,k)*xnew(OH)*xnew(NC4H10)  &
        + YCOLE(3) * rct(77,k)*xnew(OH)*xnew(C3H6) 

      L =  &
          rct(97,k)*xnew(OH) 

      xnew(ASOC_UG1E3)=  ( xold(ASOC_UG1E3) + dt2 * P) /(1.0 + dt2*L )

!-> NON_C_ASOA_NG100 

      P =   &
          0.9*rct(94,k)*xnew(OH) * xnew(ASOC_UG1 )  &
        + 1.075*rct(98,k)*xnew(OH) * xnew(NON_C_ASOA_UG1 ) 
      ! L = 0.0
 

      xnew(NON_C_ASOA_NG100)=  xold(NON_C_ASOA_NG100) + dt2 * P

!-> NON_C_ASOA_UG1 

      P =   &
          YNOXY(0) * 1.36e-11*xnew(OXYL)*xnew(OH)  &
        + YNALK(0) * rct(43,k)*xnew(OH)*xnew(NC4H10)  &
        + YNOLE(0) * rct(77,k)*xnew(OH)*xnew(C3H6)  &
        + 0.9*rct(95,k)*xnew(OH) * xnew(ASOC_UG10 )  &
        + 1.075*rct(99,k)*xnew(OH) * xnew(NON_C_ASOA_UG10 ) 

      L =  &
          rct(98,k)*xnew(OH) 

      xnew(NON_C_ASOA_UG1)=  ( xold(NON_C_ASOA_UG1) + dt2 * P) /(1.0 + dt2*L )

!-> NON_C_ASOA_UG10 

      P =   &
          YNOXY(1) * 1.36e-11*xnew(OXYL)*xnew(OH)  &
        + YNALK(1) * rct(43,k)*xnew(OH)*xnew(NC4H10)  &
        + YNOLE(1) * rct(77,k)*xnew(OH)*xnew(C3H6)  &
        + 0.9*rct(96,k)*xnew(OH) * xnew(ASOC_UG1E2 )  &
        + 1.075*rct(100,k)*xnew(OH) * xnew(NON_C_ASOA_UG1E2 ) 

      L =  &
          rct(99,k)*xnew(OH) 

      xnew(NON_C_ASOA_UG10)=  ( xold(NON_C_ASOA_UG10) + dt2 * P) /(1.0 + dt2*L )

!-> NON_C_ASOA_UG1E2 

      P =   &
          YNOXY(2) * 1.36e-11*xnew(OXYL)*xnew(OH)  &
        + YNALK(2) * rct(43,k)*xnew(OH)*xnew(NC4H10)  &
        + YNOLE(2) * rct(77,k)*xnew(OH)*xnew(C3H6)  &
        + 0.9*rct(97,k)*xnew(OH) * xnew(ASOC_UG1E3 )  &
        + 1.075*rct(101,k)*xnew(OH) * xnew(NON_C_ASOA_UG1E3 ) 

      L =  &
          rct(100,k)*xnew(OH) 

      xnew(NON_C_ASOA_UG1E2)=  ( xold(NON_C_ASOA_UG1E2) + dt2 * P) /(1.0 + dt2*L )

!-> NON_C_ASOA_UG1E3 

      P =   &
          YNOXY(3) * 1.36e-11*xnew(OXYL)*xnew(OH)  &
        + YNALK(3) * rct(43,k)*xnew(OH)*xnew(NC4H10)  &
        + YNOLE(3) * rct(77,k)*xnew(OH)*xnew(C3H6) 

      L =  &
          rct(101,k)*xnew(OH) 

      xnew(NON_C_ASOA_UG1E3)=  ( xold(NON_C_ASOA_UG1E3) + dt2 * P) /(1.0 + dt2*L )

!-> BSOC_NG100 

      P =   &
          rct(102,k)*xnew(OH) * xnew(BSOC_UG1 ) 
      ! L = 0.0
 

      xnew(BSOC_NG100)=  xold(BSOC_NG100) + dt2 * P

!-> BSOC_UG1 

      P =   &
          YCISOP(0) * rct(54,k)*xnew(OH)*xnew(C5H8)  &
        + YCTERP(0) * rct(83,k)*xnew(OH)*xnew(APINENE)  &
        + YCTERP(0) * rct(84,k)*xnew(OH)*xnew(BPINENE)  &
        + YCTERP(0) * rct(85,k)*xnew(OH)*xnew(XTERP)  &
        + YCTERP(0) * rct(88,k)*xnew(O3)*xnew(APINENE)  &
        + YCTERP(0) * rct(89,k)*xnew(O3)*xnew(BPINENE)  &
        + YCTERP(0) * rct(90,k)*xnew(O3)*xnew(XTERP)  &
        + YCTERP(0) * rct(93,k)*xnew(NO3)*xnew(APINENE)  &
        + YCTERP(0) * 2.51e-12*xnew(NO3)*xnew(BPINENE)  &
        + YCTERP(0) * rct(92,k)*xnew(NO3)*xnew(XTERP)  &
        + rct(103,k)*xnew(OH) * xnew(BSOC_UG10 ) 

      L =  &
          rct(102,k)*xnew(OH) 

      xnew(BSOC_UG1)=  ( xold(BSOC_UG1) + dt2 * P) /(1.0 + dt2*L )

!-> BSOC_UG10 

      P =   &
          YCISOP(1) * rct(54,k)*xnew(OH)*xnew(C5H8)  &
        + YCTERP(1) * rct(83,k)*xnew(OH)*xnew(APINENE)  &
        + YCTERP(1) * rct(84,k)*xnew(OH)*xnew(BPINENE)  &
        + YCTERP(1) * rct(85,k)*xnew(OH)*xnew(XTERP)  &
        + YCTERP(1) * rct(88,k)*xnew(O3)*xnew(APINENE)  &
        + YCTERP(1) * rct(89,k)*xnew(O3)*xnew(BPINENE)  &
        + YCTERP(1) * rct(90,k)*xnew(O3)*xnew(XTERP)  &
        + YCTERP(1) * rct(93,k)*xnew(NO3)*xnew(APINENE)  &
        + YCTERP(1) * 2.51e-12*xnew(NO3)*xnew(BPINENE)  &
        + YCTERP(1) * rct(92,k)*xnew(NO3)*xnew(XTERP)  &
        + rct(104,k)*xnew(OH) * xnew(BSOC_UG1E2 ) 

      L =  &
          rct(103,k)*xnew(OH) 

      xnew(BSOC_UG10)=  ( xold(BSOC_UG10) + dt2 * P) /(1.0 + dt2*L )

!-> BSOC_UG1E2 

      P =   &
          YCISOP(2) * rct(54,k)*xnew(OH)*xnew(C5H8)  &
        + YCTERP(2) * rct(83,k)*xnew(OH)*xnew(APINENE)  &
        + YCTERP(2) * rct(84,k)*xnew(OH)*xnew(BPINENE)  &
        + YCTERP(2) * rct(85,k)*xnew(OH)*xnew(XTERP)  &
        + YCTERP(2) * rct(88,k)*xnew(O3)*xnew(APINENE)  &
        + YCTERP(2) * rct(89,k)*xnew(O3)*xnew(BPINENE)  &
        + YCTERP(2) * rct(90,k)*xnew(O3)*xnew(XTERP)  &
        + YCTERP(2) * rct(93,k)*xnew(NO3)*xnew(APINENE)  &
        + YCTERP(2) * 2.51e-12*xnew(NO3)*xnew(BPINENE)  &
        + YCTERP(2) * rct(92,k)*xnew(NO3)*xnew(XTERP)  &
        + rct(105,k)*xnew(OH) * xnew(BSOC_UG1E3 ) 

      L =  &
          rct(104,k)*xnew(OH) 

      xnew(BSOC_UG1E2)=  ( xold(BSOC_UG1E2) + dt2 * P) /(1.0 + dt2*L )

!-> BSOC_UG1E3 

      P =   &
          YCISOP(3) * rct(54,k)*xnew(OH)*xnew(C5H8)  &
        + YCTERP(3) * rct(83,k)*xnew(OH)*xnew(APINENE)  &
        + YCTERP(3) * rct(84,k)*xnew(OH)*xnew(BPINENE)  &
        + YCTERP(3) * rct(85,k)*xnew(OH)*xnew(XTERP)  &
        + YCTERP(3) * rct(88,k)*xnew(O3)*xnew(APINENE)  &
        + YCTERP(3) * rct(89,k)*xnew(O3)*xnew(BPINENE)  &
        + YCTERP(3) * rct(90,k)*xnew(O3)*xnew(XTERP)  &
        + YCTERP(3) * rct(93,k)*xnew(NO3)*xnew(APINENE)  &
        + YCTERP(3) * 2.51e-12*xnew(NO3)*xnew(BPINENE)  &
        + YCTERP(3) * rct(92,k)*xnew(NO3)*xnew(XTERP) 

      L =  &
          rct(105,k)*xnew(OH) 

      xnew(BSOC_UG1E3)=  ( xold(BSOC_UG1E3) + dt2 * P) /(1.0 + dt2*L )

!-> NON_C_BSOA_NG100 

      P =   &
          0.9*rct(102,k)*xnew(OH) * xnew(BSOC_UG1 )  &
        + 1.075*rct(106,k)*xnew(OH) * xnew(NON_C_BSOA_UG1 ) 
      ! L = 0.0
 

      xnew(NON_C_BSOA_NG100)=  xold(NON_C_BSOA_NG100) + dt2 * P

!-> NON_C_BSOA_UG1 

      P =   &
          YNISOP(0) * rct(54,k)*xnew(OH)*xnew(C5H8)  &
        + YNTERP(0) * rct(83,k)*xnew(OH)*xnew(APINENE)  &
        + YNTERP(0) * rct(84,k)*xnew(OH)*xnew(BPINENE)  &
        + YNTERP(0) * rct(85,k)*xnew(OH)*xnew(XTERP)  &
        + YNTERP(0) * rct(88,k)*xnew(O3)*xnew(APINENE)  &
        + YNTERP(0) * rct(89,k)*xnew(O3)*xnew(BPINENE)  &
        + YNTERP(0) * rct(90,k)*xnew(O3)*xnew(XTERP)  &
        + YNTERP(0) * rct(93,k)*xnew(NO3)*xnew(APINENE)  &
        + YNTERP(0) * 2.51e-12*xnew(NO3)*xnew(BPINENE)  &
        + YNTERP(0) * rct(92,k)*xnew(NO3)*xnew(XTERP)  &
        + 0.9*rct(103,k)*xnew(OH) * xnew(BSOC_UG10 )  &
        + 1.075*rct(107,k)*xnew(OH) * xnew(NON_C_BSOA_UG10 ) 

      L =  &
          rct(106,k)*xnew(OH) 

      xnew(NON_C_BSOA_UG1)=  ( xold(NON_C_BSOA_UG1) + dt2 * P) /(1.0 + dt2*L )

!-> NON_C_BSOA_UG10 

      P =   &
          YNISOP(1) * rct(54,k)*xnew(OH)*xnew(C5H8)  &
        + YNTERP(1) * rct(83,k)*xnew(OH)*xnew(APINENE)  &
        + YNTERP(1) * rct(84,k)*xnew(OH)*xnew(BPINENE)  &
        + YNTERP(1) * rct(85,k)*xnew(OH)*xnew(XTERP)  &
        + YNTERP(1) * rct(88,k)*xnew(O3)*xnew(APINENE)  &
        + YNTERP(1) * rct(89,k)*xnew(O3)*xnew(BPINENE)  &
        + YNTERP(1) * rct(90,k)*xnew(O3)*xnew(XTERP)  &
        + YNTERP(1) * rct(93,k)*xnew(NO3)*xnew(APINENE)  &
        + YNTERP(1) * 2.51e-12*xnew(NO3)*xnew(BPINENE)  &
        + YNTERP(1) * rct(92,k)*xnew(NO3)*xnew(XTERP)  &
        + 0.9*rct(104,k)*xnew(OH) * xnew(BSOC_UG1E2 )  &
        + 1.075*rct(108,k)*xnew(OH) * xnew(NON_C_BSOA_UG1E2 ) 

      L =  &
          rct(107,k)*xnew(OH) 

      xnew(NON_C_BSOA_UG10)=  ( xold(NON_C_BSOA_UG10) + dt2 * P) /(1.0 + dt2*L )

!-> NON_C_BSOA_UG1E2 

      P =   &
          YNISOP(2) * rct(54,k)*xnew(OH)*xnew(C5H8)  &
        + YNTERP(2) * rct(83,k)*xnew(OH)*xnew(APINENE)  &
        + YNTERP(2) * rct(84,k)*xnew(OH)*xnew(BPINENE)  &
        + YNTERP(2) * rct(85,k)*xnew(OH)*xnew(XTERP)  &
        + YNTERP(2) * rct(88,k)*xnew(O3)*xnew(APINENE)  &
        + YNTERP(2) * rct(89,k)*xnew(O3)*xnew(BPINENE)  &
        + YNTERP(2) * rct(90,k)*xnew(O3)*xnew(XTERP)  &
        + YNTERP(2) * rct(93,k)*xnew(NO3)*xnew(APINENE)  &
        + YNTERP(2) * 2.51e-12*xnew(NO3)*xnew(BPINENE)  &
        + YNTERP(2) * rct(92,k)*xnew(NO3)*xnew(XTERP)  &
        + 0.9*rct(105,k)*xnew(OH) * xnew(BSOC_UG1E3 )  &
        + 1.075*rct(109,k)*xnew(OH) * xnew(NON_C_BSOA_UG1E3 ) 

      L =  &
          rct(108,k)*xnew(OH) 

      xnew(NON_C_BSOA_UG1E2)=  ( xold(NON_C_BSOA_UG1E2) + dt2 * P) /(1.0 + dt2*L )

!-> NON_C_BSOA_UG1E3 

      P =   &
          YNISOP(3) * rct(54,k)*xnew(OH)*xnew(C5H8)  &
        + YNTERP(3) * rct(83,k)*xnew(OH)*xnew(APINENE)  &
        + YNTERP(3) * rct(84,k)*xnew(OH)*xnew(BPINENE)  &
        + YNTERP(3) * rct(85,k)*xnew(OH)*xnew(XTERP)  &
        + YNTERP(3) * rct(88,k)*xnew(O3)*xnew(APINENE)  &
        + YNTERP(3) * rct(89,k)*xnew(O3)*xnew(BPINENE)  &
        + YNTERP(3) * rct(90,k)*xnew(O3)*xnew(XTERP)  &
        + YNTERP(3) * rct(93,k)*xnew(NO3)*xnew(APINENE)  &
        + YNTERP(3) * 2.51e-12*xnew(NO3)*xnew(BPINENE)  &
        + YNTERP(3) * rct(92,k)*xnew(NO3)*xnew(XTERP) 

      L =  &
          rct(109,k)*xnew(OH) 

      xnew(NON_C_BSOA_UG1E3)=  ( xold(NON_C_BSOA_UG1E3) + dt2 * P) /(1.0 + dt2*L )

!-> FFFUEL_NG10 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for FFFUEL_NG10! xnew(FFFUEL_NG10)=  max(0.0, xold(FFFUEL_NG10))

!-> WOODOA_NG10 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for WOODOA_NG10! xnew(WOODOA_NG10)=  max(0.0, xold(WOODOA_NG10))

!-> FFIREOA_NG10 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for FFIREOA_NG10! xnew(FFIREOA_NG10)=  max(0.0, xold(FFIREOA_NG10))

!-> SEASALT_F 

      P =   &
          rcemis(SEASALT_F,k) 
      ! L = 0.0
 

      xnew(SEASALT_F)=  xold(SEASALT_F) + dt2 * P

!-> SEASALT_C 

      P =   &
          rcemis(SEASALT_C,k) 
      ! L = 0.0
 

      xnew(SEASALT_C)=  xold(SEASALT_C) + dt2 * P

!-> DUST_ROAD_F 

      P =   &
          rcemis(DUST_ROAD_F,k) 
      ! L = 0.0
 

      xnew(DUST_ROAD_F)=  xold(DUST_ROAD_F) + dt2 * P

!-> DUST_ROAD_C 

      P =   &
          rcemis(DUST_ROAD_C,k) 
      ! L = 0.0
 

      xnew(DUST_ROAD_C)=  xold(DUST_ROAD_C) + dt2 * P

!-> DUST_WB_F 

      P =   &
          rcemis(DUST_WB_F,k) 
      ! L = 0.0
 

      xnew(DUST_WB_F)=  xold(DUST_WB_F) + dt2 * P

!-> DUST_WB_C 

      P =   &
          rcemis(DUST_WB_C,k) 
      ! L = 0.0
 

      xnew(DUST_WB_C)=  xold(DUST_WB_C) + dt2 * P

!-> DUST_SAH_F 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for DUST_SAH_F! xnew(DUST_SAH_F)=  max(0.0, xold(DUST_SAH_F))

!-> DUST_SAH_C 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for DUST_SAH_C! xnew(DUST_SAH_C)=  max(0.0, xold(DUST_SAH_C))

!-> RN222 

      P =   &
          rcemis(RN222,k) 

      L =  &
          2.1e-6 

      xnew(RN222)=  ( xold(RN222) + dt2 * P) /(1.0 + dt2*L )

!-> RNWATER 
      ! P = 0.0
 

      L =  &
          2.1e-6 

      xnew(RNWATER)=  xold(RNWATER) / ( 1.0 + dt2 * L )

!-> PB210 

      P =   &
          2.1e-6 * xnew(RN222 ) 
      ! L = 0.0
 

      xnew(PB210)=  xold(PB210) + dt2 * P
