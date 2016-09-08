
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
          rcemis(NH3,k) 
      ! L = 0.0
 

      xnew(NH3)=  xold(NH3) + dt2 * P

!-> NO3_F 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for NO3_F! xnew(NO3_F)=  max(0.0, xold(NO3_F))

!-> NO3_C 

      P =   &
          rct(71,k) * xnew(HNO3 )  &
        + rct(72,k) * xnew(HNO3 ) 
      ! L = 0.0
 

      xnew(NO3_C)=  xold(NO3_C) + dt2 * P

!-> NH4_F 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for NH4_F! xnew(NH4_F)=  max(0.0, xold(NH4_F))

!-> DUMMY 
      ! P = 0.0
 

      L =  &
          rct(75,k) 

      xnew(DUMMY)=  xold(DUMMY) / ( 1.0 + dt2 * L )

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
          rct(96,k) 

      xnew(EC_F_WOOD_NEW)=  ( xold(EC_F_WOOD_NEW) + dt2 * P) /(1.0 + dt2*L )

!-> EC_F_WOOD_AGE 

      P =   &
          rcemis(EC_F_WOOD_AGE,k)  &
        + rct(96,k) * xnew(EC_F_WOOD_NEW ) 
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
          rct(97,k) 

      xnew(EC_F_FFUEL_NEW)=  ( xold(EC_F_FFUEL_NEW) + dt2 * P) /(1.0 + dt2*L )

!-> EC_F_FFUEL_AGE 

      P =   &
          rcemis(EC_F_FFUEL_AGE,k)  &
        + rct(97,k) * xnew(EC_F_FFUEL_NEW ) 
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

!-> ASOC_NG100 

      P =   &
          rct(80,k)*xnew(OH) * xnew(ASOC_UG1 ) 
      ! L = 0.0
 

      xnew(ASOC_NG100)=  xold(ASOC_NG100) + dt2 * P

!-> ASOC_UG1 

      P =   &
          0.00206*rct(40,k)*xnew(NO)*xnew(PRRO2)  &
        + 0.010294*rct(45,k)*xnew(PRRO2)*xnew(HO2)  &
        + 0.008413*rct(40,k)*xnew(OXYO2)*xnew(NO)  &
        + 0.315476*rct(46,k)*xnew(OXYO2)*xnew(HO2)  &
        + rct(81,k)*xnew(OH) * xnew(ASOC_UG10 ) 

      L =  &
          rct(80,k)*xnew(OH) 

      xnew(ASOC_UG1)=  ( xold(ASOC_UG1) + dt2 * P) /(1.0 + dt2*L )

!-> ASOC_UG10 

      P =   &
          0.108*rct(40,k)*xnew(NO)*xnew(SECC4H9O2)  &
        + 0.2132*rct(41,k)*xnew(SECC4H9O2)*xnew(HO2)  &
        + 0.01029*rct(40,k)*xnew(NO)*xnew(PRRO2)  &
        + 0.018529*rct(45,k)*xnew(PRRO2)*xnew(HO2)  &
        + 0.820238*rct(40,k)*xnew(OXYO2)*xnew(NO)  &
        + 1.261905*rct(46,k)*xnew(OXYO2)*xnew(HO2)  &
        + rct(82,k)*xnew(OH) * xnew(ASOC_UG1E2 ) 

      L =  &
          rct(81,k)*xnew(OH) 

      xnew(ASOC_UG10)=  ( xold(ASOC_UG10) + dt2 * P) /(1.0 + dt2*L )

!-> ASOC_UG1E2 

      P =   &
          0.078235*rct(40,k)*xnew(NO)*xnew(PRRO2)  &
        + 0.123529*rct(45,k)*xnew(PRRO2)*xnew(HO2)  &
        + 1.261905*rct(40,k)*xnew(OXYO2)*xnew(NO)  &
        + 1.577381*rct(46,k)*xnew(OXYO2)*xnew(HO2)  &
        + rct(83,k)*xnew(OH) * xnew(ASOC_UG1E3 ) 

      L =  &
          rct(82,k)*xnew(OH) 

      xnew(ASOC_UG1E2)=  ( xold(ASOC_UG1E2) + dt2 * P) /(1.0 + dt2*L )

!-> ASOC_UG1E3 

      P =   &
          0.30882*rct(40,k)*xnew(NO)*xnew(PRRO2)  &
        + 0.463235*rct(45,k)*xnew(PRRO2)*xnew(HO2)  &
        + 1.829762*rct(40,k)*xnew(OXYO2)*xnew(NO)  &
        + 2.208333*rct(46,k)*xnew(OXYO2)*xnew(HO2) 

      L =  &
          rct(83,k)*xnew(OH) 

      xnew(ASOC_UG1E3)=  ( xold(ASOC_UG1E3) + dt2 * P) /(1.0 + dt2*L )

!-> NON_C_ASOA_NG100 

      P =   &
          0.9*rct(80,k)*xnew(OH) * xnew(ASOC_UG1 )  &
        + 1.075*rct(84,k)*xnew(OH) * xnew(NON_C_ASOA_UG1 ) 
      ! L = 0.0
 

      xnew(NON_C_ASOA_NG100)=  xold(NON_C_ASOA_NG100) + dt2 * P

!-> NON_C_ASOA_UG1 

      P =   &
          0.0173*rct(40,k)*xnew(NO)*xnew(PRRO2)  &
        + 0.08647*rct(45,k)*xnew(PRRO2)*xnew(HO2)  &
        + 0.11105*rct(40,k)*xnew(OXYO2)*xnew(NO)  &
        + 4.16429*rct(46,k)*xnew(OXYO2)*xnew(HO2)  &
        + 0.9*rct(81,k)*xnew(OH) * xnew(ASOC_UG10 )  &
        + 1.075*rct(85,k)*xnew(OH) * xnew(NON_C_ASOA_UG10 ) 

      L =  &
          rct(84,k)*xnew(OH) 

      xnew(NON_C_ASOA_UG1)=  ( xold(NON_C_ASOA_UG1) + dt2 * P) /(1.0 + dt2*L )

!-> NON_C_ASOA_UG10 

      P =   &
          0.90753*rct(40,k)*xnew(NO)*xnew(SECC4H9O2)  &
        + 1.7912*rct(41,k)*xnew(SECC4H9O2)*xnew(HO2)  &
        + 0.08647*rct(40,k)*xnew(NO)*xnew(PRRO2)  &
        + 0.15565*rct(45,k)*xnew(PRRO2)*xnew(HO2)  &
        + 10.8271*rct(40,k)*xnew(OXYO2)*xnew(NO)  &
        + 16.6571*rct(46,k)*xnew(OXYO2)*xnew(HO2)  &
        + 0.9*rct(82,k)*xnew(OH) * xnew(ASOC_UG1E2 )  &
        + 1.075*rct(86,k)*xnew(OH) * xnew(NON_C_ASOA_UG1E2 ) 

      L =  &
          rct(85,k)*xnew(OH) 

      xnew(NON_C_ASOA_UG10)=  ( xold(NON_C_ASOA_UG10) + dt2 * P) /(1.0 + dt2*L )

!-> NON_C_ASOA_UG1E2 

      P =   &
          0.65718*rct(40,k)*xnew(NO)*xnew(PRRO2)  &
        + 1.03765*rct(45,k)*xnew(PRRO2)*xnew(HO2)  &
        + 16.6571*rct(40,k)*xnew(OXYO2)*xnew(NO)  &
        + 20.82143*rct(46,k)*xnew(OXYO2)*xnew(HO2)  &
        + 0.9*rct(83,k)*xnew(OH) * xnew(ASOC_UG1E3 )  &
        + 1.075*rct(87,k)*xnew(OH) * xnew(NON_C_ASOA_UG1E3 ) 

      L =  &
          rct(86,k)*xnew(OH) 

      xnew(NON_C_ASOA_UG1E2)=  ( xold(NON_C_ASOA_UG1E2) + dt2 * P) /(1.0 + dt2*L )

!-> NON_C_ASOA_UG1E3 

      P =   &
          2.59412*rct(40,k)*xnew(NO)*xnew(PRRO2)  &
        + 3.89118*rct(45,k)*xnew(PRRO2)*xnew(HO2)  &
        + 24.15286*rct(40,k)*xnew(OXYO2)*xnew(NO)  &
        + 29.15*rct(46,k)*xnew(OXYO2)*xnew(HO2) 

      L =  &
          rct(87,k)*xnew(OH) 

      xnew(NON_C_ASOA_UG1E3)=  ( xold(NON_C_ASOA_UG1E3) + dt2 * P) /(1.0 + dt2*L )

!-> BSOC_NG100 

      P =   &
          rct(88,k)*xnew(OH) * xnew(BSOC_UG1 ) 
      ! L = 0.0
 

      xnew(BSOC_NG100)=  xold(BSOC_NG100) + dt2 * P

!-> BSOC_UG1 

      P =   &
          0.002833*rct(40,k)*xnew(ISRO2)*xnew(NO)  &
        + 0.0255*rct(47,k)*xnew(ISRO2)*xnew(HO2)  &
        + 0.08*rct(40,k)*xnew(NO) * xnew(TERPPEROXY )  &
        + 0.715333*rct(79,k)*xnew(HO2) * xnew(TERPPEROXY )  &
        + rct(89,k)*xnew(OH) * xnew(BSOC_UG10 ) 

      L =  &
          rct(88,k)*xnew(OH) 

      xnew(BSOC_UG1)=  ( xold(BSOC_UG1) + dt2 * P) /(1.0 + dt2*L )

!-> BSOC_UG10 

      P =   &
          0.065167*rct(40,k)*xnew(ISRO2)*xnew(NO)  &
        + 0.085*rct(47,k)*xnew(ISRO2)*xnew(HO2)  &
        + 0.813333*rct(40,k)*xnew(NO) * xnew(TERPPEROXY )  &
        + 0.612*rct(79,k)*xnew(HO2) * xnew(TERPPEROXY )  &
        + rct(90,k)*xnew(OH) * xnew(BSOC_UG1E2 ) 

      L =  &
          rct(89,k)*xnew(OH) 

      xnew(BSOC_UG10)=  ( xold(BSOC_UG10) + dt2 * P) /(1.0 + dt2*L )

!-> BSOC_UG1E2 

      P =   &
          0.0425*rct(40,k)*xnew(ISRO2)*xnew(NO)  &
        + 0.0425*rct(47,k)*xnew(ISRO2)*xnew(HO2)  &
        + 1.34*rct(40,k)*xnew(NO) * xnew(TERPPEROXY )  &
        + 2.391333*rct(79,k)*xnew(HO2) * xnew(TERPPEROXY )  &
        + rct(91,k)*xnew(OH) * xnew(BSOC_UG1E3 ) 

      L =  &
          rct(90,k)*xnew(OH) 

      xnew(BSOC_UG1E2)=  ( xold(BSOC_UG1E2) + dt2 * P) /(1.0 + dt2*L )

!-> BSOC_UG1E3 

      P =   &
          3.333333*rct(40,k)*xnew(NO) * xnew(TERPPEROXY )  &
        + 4.05*rct(79,k)*xnew(HO2) * xnew(TERPPEROXY ) 

      L =  &
          rct(91,k)*xnew(OH) 

      xnew(BSOC_UG1E3)=  ( xold(BSOC_UG1E3) + dt2 * P) /(1.0 + dt2*L )

!-> NON_C_BSOA_NG100 

      P =   &
          0.9*rct(88,k)*xnew(OH) * xnew(BSOC_UG1 )  &
        + 1.075*rct(92,k)*xnew(OH) * xnew(NON_C_BSOA_UG1 ) 
      ! L = 0.0
 

      xnew(NON_C_BSOA_NG100)=  xold(NON_C_BSOA_NG100) + dt2 * P

!-> NON_C_BSOA_UG1 

      P =   &
          0.034*rct(40,k)*xnew(ISRO2)*xnew(NO)  &
        + 0.306*rct(47,k)*xnew(ISRO2)*xnew(HO2)  &
        + 0.672*rct(40,k)*xnew(NO) * xnew(TERPPEROXY )  &
        + 6.009*rct(79,k)*xnew(HO2) * xnew(TERPPEROXY )  &
        + 0.9*rct(89,k)*xnew(OH) * xnew(BSOC_UG10 )  &
        + 1.075*rct(93,k)*xnew(OH) * xnew(NON_C_BSOA_UG10 ) 

      L =  &
          rct(92,k)*xnew(OH) 

      xnew(NON_C_BSOA_UG1)=  ( xold(NON_C_BSOA_UG1) + dt2 * P) /(1.0 + dt2*L )

!-> NON_C_BSOA_UG10 

      P =   &
          0.782*rct(40,k)*xnew(ISRO2)*xnew(NO)  &
        + 1.02*rct(47,k)*xnew(ISRO2)*xnew(HO2)  &
        + 6.832*rct(40,k)*xnew(NO) * xnew(TERPPEROXY )  &
        + 5.1408*rct(79,k)*xnew(HO2) * xnew(TERPPEROXY )  &
        + 0.9*rct(90,k)*xnew(OH) * xnew(BSOC_UG1E2 )  &
        + 1.075*rct(94,k)*xnew(OH) * xnew(NON_C_BSOA_UG1E2 ) 

      L =  &
          rct(93,k)*xnew(OH) 

      xnew(NON_C_BSOA_UG10)=  ( xold(NON_C_BSOA_UG10) + dt2 * P) /(1.0 + dt2*L )

!-> NON_C_BSOA_UG1E2 

      P =   &
          0.51*rct(40,k)*xnew(ISRO2)*xnew(NO)  &
        + 0.51*rct(47,k)*xnew(ISRO2)*xnew(HO2)  &
        + 11.256*rct(40,k)*xnew(NO) * xnew(TERPPEROXY )  &
        + 20.0872*rct(79,k)*xnew(HO2) * xnew(TERPPEROXY )  &
        + 0.9*rct(91,k)*xnew(OH) * xnew(BSOC_UG1E3 )  &
        + 1.075*rct(95,k)*xnew(OH) * xnew(NON_C_BSOA_UG1E3 ) 

      L =  &
          rct(94,k)*xnew(OH) 

      xnew(NON_C_BSOA_UG1E2)=  ( xold(NON_C_BSOA_UG1E2) + dt2 * P) /(1.0 + dt2*L )

!-> NON_C_BSOA_UG1E3 

      P =   &
          28.*rct(40,k)*xnew(NO) * xnew(TERPPEROXY )  &
        + 34.02*rct(79,k)*xnew(HO2) * xnew(TERPPEROXY ) 

      L =  &
          rct(95,k)*xnew(OH) 

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
