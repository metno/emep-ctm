
!-> OD 

      P =   &
          rcphot(IDBO3,K) * xnew(O3 ) 

      L =  &
          rct(2,k)  &
        + rct(3,k)  &
        + rct(4,k) 

      xnew(OD)=  ( xold(OD) + dt2 * P) /(1.0 + dt2*L )

!-> OP 

      P =   &
          rct(2,k) * xnew(OD )  &
        + rct(3,k) * xnew(OD )  &
        + 0.3*rct(50,k) * xnew(C5H8 ) * xnew(O3 )  &
        + 0.2*rct(58,k) * xnew(MVK ) * xnew(O3 )  &
        + rcphot(IDAO3,K) * xnew(O3 )  &
        + rcphot(IDNO2,K) * xnew(NO2 )  &
        + rcphot(IDNO3,K) * xnew(NO3 ) 

      L =  &
          rct(1,k)  &
        + rct(63,k)* xnew(NO ) 

      xnew(OP)=  ( xold(OP) + dt2 * P) /(1.0 + dt2*L )

!-> OH 

      P =   &
          2.*rct(4,k) * xnew(OD )  &
        + rct(8,k) * xnew(O3 ) * xnew(HO2 )  &
        + rct(10,k) * xnew(NO ) * xnew(HO2 )  &
        + 0.4*rct(26,k) * xnew(CH3O2H ) * xnew(OH )  &
        + 8.01e-12 * xnew(C2H5OOH ) * xnew(OH )  &
        + 0.44*rct(37,k) * xnew(CH3COO2 ) * xnew(HO2 )  &
        + 2.15e-11 * xnew(BURO2H ) * xnew(OH )  &
        + 1.38e-11 * xnew(ETRO2H ) * xnew(OH )  &
        + 0.13*rct(43,k) * xnew(C2H4 ) * xnew(O3 )  &
        + 0.36*rct(44,k) * xnew(O3 ) * xnew(C3H6 )  &
        + 2.44e-11 * xnew(PRRO2H ) * xnew(OH )  &
        + 0.55*rct(50,k) * xnew(C5H8 ) * xnew(O3 )  &
        + 7.5e-11 * xnew(ISRO2H ) * xnew(OH )  &
        + 0.82*rct(53,k) * xnew(MACR ) * xnew(O3 )  &
        + 0.08*rct(58,k) * xnew(MVK ) * xnew(O3 )  &
        + 2*rcphot(IDH2O2,K) * xnew(H2O2 )  &
        + rcphot(IDHNO3,K) * xnew(HNO3 )  &
        + rcphot(IDCH3O2H,K) * xnew(CH3O2H )  &
        + rcphot(IDCH3O2H,K) * xnew(C2H5OOH )  &
        + rcphot(IDCH3O2H,K) * xnew(ETRO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(BURO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(MEKO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(CH3COO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(OXYO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(MALO2H )  &
        + 0.222*rcphot(IDNO2,K) * xnew(HONO ) 

      L =  &
          rct(7,k)* xnew(O3 )  &
        + rct(12,k)* xnew(HO2 )  &
        + rct(13,k)* xnew(H2O2 )  &
        + rct(14,k)* xnew(H2 )  &
        + rct(15,k)* xnew(HNO3 )  &
        + rct(18,k)* xnew(HONO )  &
        + 2e-12*AQRCK(ICLOHSO2,K)* xnew(SO2 )  &
        + rct(19,k)* xnew(CH4 )  &
        + rct(20,k)* xnew(CO )  &
        + rct(24,k)* xnew(CH3OH )  &
        + rct(26,k)* xnew(CH3O2H )  &
        + rct(27,k)* xnew(HCHO )  &
        + rct(29,k)* xnew(C2H6 )  &
        + 8.01e-12* xnew(C2H5OOH )  &
        + rct(32,k)* xnew(C2H5OOH )  &
        + rct(33,k)* xnew(CH3CHO )  &
        + rct(32,k)* xnew(CH3COO2H )  &
        + rct(38,k)* xnew(C2H5OH )  &
        + rct(39,k)* xnew(NC4H10 )  &
        + rct(42,k)* xnew(MEK )  &
        + rct(32,k)* xnew(MEKO2H )  &
        + rct(32,k)* xnew(BURO2H )  &
        + 2.15e-11* xnew(BURO2H )  &
        + 1.38e-11* xnew(ETRO2H )  &
        + rct(32,k)* xnew(ETRO2H )  &
        + 2.44e-11* xnew(PRRO2H )  &
        + rct(32,k)* xnew(PRRO2H )  &
        + 1.36e-11* xnew(OXYL )  &
        + 4.2e-11* xnew(OXYO2H )  &
        + 5.58e-11* xnew(MAL )  &
        + rct(32,k)* xnew(MALO2H )  &
        + rct(48,k)* xnew(GLYOX )  &
        + rct(49,k)* xnew(MGLYOX )  &
        + rct(51,k)* xnew(C5H8 )  &
        + rct(52,k)* xnew(MVK )  &
        + 7.5e-11* xnew(ISRO2H )  &
        + rct(54,k)* xnew(MACR )  &
        + 2.82e-11* xnew(MACROOH )  &
        + rct(56,k)* xnew(ACETOL )  &
        + 5.96e-11* xnew(ISNI )  &
        + 1.87e-11* xnew(MACO3H )  &
        + 1.51e-11* xnew(MACO2H )  &
        + 2.0e-11* xnew(ISONO3H )  &
        + 2.2e-11* xnew(MVKO2H )  &
        + 3.7e-11* xnew(ISNIRH )  &
        + 2.9e-11* xnew(MPAN )  &
        + rct(66,k)* xnew(NO2 )  &
        + rct(69,k)* xnew(C2H4 )  &
        + rct(70,k)* xnew(C3H6 )  &
        + rct(71,k)* xnew(NO ) 

      xnew(OH)=  ( xold(OH) + dt2 * P) /(1.0 + dt2*L )

!-> HO2 

      P =   &
          rct(7,k) * xnew(O3 ) * xnew(OH )  &
        + rct(13,k) * xnew(OH ) * xnew(H2O2 )  &
        + rct(14,k) * xnew(OH ) * xnew(H2 )  &
        + 2e-12*AQRCK(ICLOHSO2,K) * xnew(OH ) * xnew(SO2 )  &
        + rct(20,k) * xnew(OH ) * xnew(CO )  &
        + rct(21,k) * xnew(CH3O2 ) * xnew(NO )  &
        + 2.*rct(22,k) * xnew(CH3O2 ) * xnew(CH3O2 )  &
        + rct(24,k) * xnew(OH ) * xnew(CH3OH )  &
        + rct(27,k) * xnew(OH ) * xnew(HCHO )  &
        + rct(28,k) * xnew(NO3 ) * xnew(HCHO )  &
        + rct(30,k) * xnew(C2H5O2 ) * xnew(NO )  &
        + 0.9*rct(35,k) * xnew(CH3O2 ) * xnew(CH3COO2 )  &
        + rct(38,k) * xnew(OH ) * xnew(C2H5OH )  &
        + 0.65*rct(40,k) * xnew(NO ) * xnew(SECC4H9O2 )  &
        + rct(40,k) * xnew(ETRO2 ) * xnew(NO )  &
        + 0.13*rct(43,k) * xnew(C2H4 ) * xnew(O3 )  &
        + 0.28*rct(44,k) * xnew(O3 ) * xnew(C3H6 )  &
        + rct(40,k) * xnew(NO ) * xnew(PRRO2 )  &
        + rct(40,k) * xnew(OXYO2 ) * xnew(NO )  &
        + rct(40,k) * xnew(MALO2 ) * xnew(NO )  &
        + rct(48,k) * xnew(OH ) * xnew(GLYOX )  &
        + 0.06*rct(50,k) * xnew(C5H8 ) * xnew(O3 )  &
        + 0.78*rct(40,k) * xnew(ISRO2 ) * xnew(NO )  &
        + 0.95*rct(40,k) * xnew(MVKO2 ) * xnew(NO )  &
        + 0.41*rct(53,k) * xnew(MACR ) * xnew(O3 )  &
        + 0.95*rct(55,k) * xnew(MACRO2 ) * xnew(NO )  &
        + 2.5e-12 * xnew(MACRO2 ) * xnew(NO3 )  &
        + rct(56,k) * xnew(ACETOL ) * xnew(OH )  &
        + 0.06*rct(58,k) * xnew(MVK ) * xnew(O3 )  &
        + 0.05*rct(40,k) * xnew(ISNIR ) * xnew(NO )  &
        + 0.8*rct(40,k) * xnew(ISONO3 ) * xnew(NO )  &
        + 2*rcphot(IDACH2O,K) * xnew(HCHO )  &
        + rcphot(IDCH3CHO,K) * xnew(CH3CHO )  &
        + rcphot(IDCH3O2H,K) * xnew(CH3O2H )  &
        + 0.5*rcphot(IDHCOHCO,K) * xnew(GLYOX )  &
        + rcphot(IDRCOHCO,K) * xnew(MGLYOX )  &
        + rcphot(IDCH3O2H,K) * xnew(C2H5OOH )  &
        + rcphot(IDCH3O2H,K) * xnew(ETRO2H )  &
        + 0.65*rcphot(IDCH3O2H,K) * xnew(BURO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(PRRO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(OXYO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(MALO2H ) 

      L =  &
          rct(8,k)* xnew(O3 )  &
        + rct(10,k)* xnew(NO )  &
        + rct(12,k)* xnew(OH )  &
        + rct(16,k)* xnew(HO2 )  &
        + rct(16,k)* xnew(HO2 )  &
        + rct(17,k)* xnew(HO2 )  &
        + rct(17,k)* xnew(HO2 )  &
        + rct(25,k)* xnew(CH3O2 )  &
        + rct(31,k)* xnew(C2H5O2 )  &
        + rct(37,k)* xnew(CH3COO2 )  &
        + rct(41,k)* xnew(SECC4H9O2 )  &
        + rct(41,k)* xnew(MEKO2 )  &
        + 1.2e-11* xnew(ETRO2 )  &
        + rct(45,k)* xnew(PRRO2 )  &
        + rct(46,k)* xnew(OXYO2 )  &
        + rct(47,k)* xnew(MALO2 )  &
        + rct(47,k)* xnew(ISRO2 )  &
        + rct(41,k)* xnew(MACRO2 )  &
        + rct(47,k)* xnew(ISONO3 )  &
        + rct(41,k)* xnew(MVKO2 )  &
        + rct(60,k)* xnew(MACO3 )  &
        + rct(47,k)* xnew(ISNIR ) 

      xnew(HO2)=  ( xold(HO2) + dt2 * P) /(1.0 + dt2*L )

!-> CH3O2 

      P =   &
          rct(19,k) * xnew(OH ) * xnew(CH4 )  &
        + 0.6*rct(26,k) * xnew(CH3O2H ) * xnew(OH )  &
        + 0.05*rct(33,k) * xnew(OH ) * xnew(CH3CHO )  &
        + rct(34,k) * xnew(CH3COO2 ) * xnew(NO )  &
        + 0.9*rct(35,k) * xnew(CH3O2 ) * xnew(CH3COO2 )  &
        + rct(36,k) * xnew(CH3COO2 ) * xnew(CH3COO2 )  &
        + rct(36,k) * xnew(CH3COO2 ) * xnew(CH3COO2 )  &
        + 0.44*rct(37,k) * xnew(CH3COO2 ) * xnew(HO2 )  &
        + 0.28*rct(44,k) * xnew(O3 ) * xnew(C3H6 )  &
        + rcphot(IDCH3CHO,K) * xnew(CH3CHO )  &
        + rcphot(IDCH3O2H,K) * xnew(CH3COO2H ) 

      L =  &
          rct(21,k)* xnew(NO )  &
        + rct(22,k)* xnew(CH3O2 )  &
        + rct(22,k)* xnew(CH3O2 )  &
        + rct(23,k)* xnew(CH3O2 )  &
        + rct(23,k)* xnew(CH3O2 )  &
        + rct(25,k)* xnew(HO2 )  &
        + rct(35,k)* xnew(CH3COO2 ) 

      xnew(CH3O2)=  ( xold(CH3O2) + dt2 * P) /(1.0 + dt2*L )

!-> C2H5O2 

      P =   &
          rct(29,k) * xnew(OH ) * xnew(C2H6 )  &
        + rct(32,k) * xnew(C2H5OOH ) * xnew(OH )  &
        + 0.35*rct(40,k) * xnew(NO ) * xnew(SECC4H9O2 )  &
        + 0.35*rcphot(IDCH3O2H,K) * xnew(BURO2H )  &
        + rcphot(IDCH3COX,K) * xnew(MEK ) 

      L =  &
          rct(30,k)* xnew(NO )  &
        + rct(31,k)* xnew(HO2 ) 

      xnew(C2H5O2)=  ( xold(C2H5O2) + dt2 * P) /(1.0 + dt2*L )

!-> SECC4H9O2 

      P =   &
          rct(39,k) * xnew(OH ) * xnew(NC4H10 )  &
        + rct(32,k) * xnew(BURO2H ) * xnew(OH ) 

      L =  &
          rct(40,k)* xnew(NO )  &
        + rct(41,k)* xnew(HO2 ) 

      xnew(SECC4H9O2)=  ( xold(SECC4H9O2) + dt2 * P) /(1.0 + dt2*L )

!-> ISRO2 

      P =   &
          rct(51,k) * xnew(C5H8 ) * xnew(OH )  &
        + 0.12*rct(40,k) * xnew(ISRO2 ) * xnew(NO )  &
        + 7.5e-11 * xnew(ISRO2H ) * xnew(OH ) 

      L =  &
          rct(40,k)* xnew(NO )  &
        + rct(47,k)* xnew(HO2 ) 

      xnew(ISRO2)=  ( xold(ISRO2) + dt2 * P) /(1.0 + dt2*L )

!-> ETRO2 

      P =   &
          rct(32,k) * xnew(ETRO2H ) * xnew(OH )  &
        + rct(69,k) * xnew(C2H4 ) * xnew(OH ) 

      L =  &
          rct(40,k)* xnew(NO )  &
        + 1.2e-11* xnew(HO2 ) 

      xnew(ETRO2)=  ( xold(ETRO2) + dt2 * P) /(1.0 + dt2*L )

!-> PRRO2 

      P =   &
          rct(32,k) * xnew(PRRO2H ) * xnew(OH )  &
        + rct(70,k) * xnew(OH ) * xnew(C3H6 ) 

      L =  &
          rct(40,k)* xnew(NO )  &
        + rct(45,k)* xnew(HO2 ) 

      xnew(PRRO2)=  ( xold(PRRO2) + dt2 * P) /(1.0 + dt2*L )

!-> OXYO2 

      P =   &
          1.36e-11 * xnew(OXYL ) * xnew(OH )  &
        + 4.2e-11 * xnew(OXYO2H ) * xnew(OH ) 

      L =  &
          rct(40,k)* xnew(NO )  &
        + rct(46,k)* xnew(HO2 ) 

      xnew(OXYO2)=  ( xold(OXYO2) + dt2 * P) /(1.0 + dt2*L )

!-> MEKO2 

      P =   &
          rct(42,k) * xnew(OH ) * xnew(MEK )  &
        + rct(32,k) * xnew(MEKO2H ) * xnew(OH ) 

      L =  &
          rct(40,k)* xnew(NO )  &
        + rct(41,k)* xnew(HO2 ) 

      xnew(MEKO2)=  ( xold(MEKO2) + dt2 * P) /(1.0 + dt2*L )

!-> MALO2 

      P =   &
          5.58e-11 * xnew(MAL ) * xnew(OH )  &
        + rct(32,k) * xnew(MALO2H ) * xnew(OH ) 

      L =  &
          rct(40,k)* xnew(NO )  &
        + rct(47,k)* xnew(HO2 ) 

      xnew(MALO2)=  ( xold(MALO2) + dt2 * P) /(1.0 + dt2*L )

!-> MVKO2 

      P =   &
          rct(52,k) * xnew(MVK ) * xnew(OH )  &
        + 2.2e-11 * xnew(MVKO2H ) * xnew(OH ) 

      L =  &
          rct(40,k)* xnew(NO )  &
        + rct(41,k)* xnew(HO2 ) 

      xnew(MVKO2)=  ( xold(MVKO2) + dt2 * P) /(1.0 + dt2*L )

!-> MACRO2 

      P =   &
          0.5*rct(54,k) * xnew(MACR ) * xnew(OH )  &
        + 2.82e-11 * xnew(MACROOH ) * xnew(OH ) 

      L =  &
          rct(55,k)* xnew(NO )  &
        + 2.5e-12* xnew(NO3 )  &
        + rct(41,k)* xnew(HO2 ) 

      xnew(MACRO2)=  ( xold(MACRO2) + dt2 * P) /(1.0 + dt2*L )

!-> MACO3 

      P =   &
          0.5*rct(54,k) * xnew(MACR ) * xnew(OH )  &
        + 1.87e-11 * xnew(MACO3H ) * xnew(OH )  &
        + rct(68,k) * xnew(MPAN ) 

      L =  &
          rct(57,k)* xnew(NO )  &
        + rct(60,k)* xnew(HO2 )  &
        + rct(67,k)* xnew(NO2 ) 

      xnew(MACO3)=  ( xold(MACO3) + dt2 * P) /(1.0 + dt2*L )

!-> TERPPEROXY 

      P =   &
          rct(72,k)*xnew(O3) * xnew(APINENE )  &
        + rct(73,k)*xnew(OH) * xnew(APINENE )  &
        + rct(74,k)*xnew(NO3) * xnew(APINENE ) 

      L =  &
          rct(40,k)*xnew(NO)  &
        + rct(75,k)*xnew(HO2) 

      xnew(TERPPEROXY)=  ( xold(TERPPEROXY) + dt2 * P) /(1.0 + dt2*L )

!-> O3 

      P =   &
          rct(1,k) * xnew(OP )  &
        + 0.15*rct(37,k) * xnew(CH3COO2 ) * xnew(HO2 )  &
        + 0.29*rct(60,k) * xnew(MACO3 ) * xnew(HO2 ) 

      L =  &
          rct(5,k)* xnew(NO )  &
        + rct(6,k)* xnew(NO2 )  &
        + rct(7,k)* xnew(OH )  &
        + rct(8,k)* xnew(HO2 )  &
        + AQRCK(ICLRC2,K)* xnew(SO2 )  &
        + rct(43,k)* xnew(C2H4 )  &
        + rct(44,k)* xnew(C3H6 )  &
        + rct(50,k)* xnew(C5H8 )  &
        + rct(53,k)* xnew(MACR )  &
        + rct(58,k)* xnew(MVK )  &
        + rcphot(IDAO3,K)  &
        + rcphot(IDBO3,K) 

      xnew(O3)=  ( xold(O3) + dt2 * P) /(1.0 + dt2*L )

!-> NO 

      P =   &
          rct(11,k) * xnew(NO2 ) * xnew(NO3 )  &
        + rcphot(IDNO2,K) * xnew(NO2 )  &
        + 0.222*rcphot(IDNO2,K) * xnew(HONO )  &
        + rcemis(NO,k)  &
        + 0 !Skip bio rate since rcemis exists 

      L =  &
          rct(5,k)* xnew(O3 )  &
        + rct(9,k)* xnew(NO3 )  &
        + rct(10,k)* xnew(HO2 )  &
        + rct(21,k)* xnew(CH3O2 )  &
        + rct(30,k)* xnew(C2H5O2 )  &
        + rct(34,k)* xnew(CH3COO2 )  &
        + rct(40,k)* xnew(SECC4H9O2 )  &
        + rct(40,k)* xnew(MEKO2 )  &
        + rct(40,k)* xnew(ETRO2 )  &
        + rct(40,k)* xnew(PRRO2 )  &
        + rct(40,k)* xnew(OXYO2 )  &
        + rct(40,k)* xnew(MALO2 )  &
        + rct(40,k)* xnew(ISRO2 )  &
        + rct(40,k)* xnew(MVKO2 )  &
        + rct(55,k)* xnew(MACRO2 )  &
        + rct(57,k)* xnew(MACO3 )  &
        + rct(40,k)* xnew(ISNIR )  &
        + rct(40,k)* xnew(ISONO3 )  &
        + rct(63,k)* xnew(OP )  &
        + rct(71,k)* xnew(OH ) 

      xnew(NO)=  ( xold(NO) + dt2 * P) /(1.0 + dt2*L )

!-> NO2 

      P =   &
          rct(5,k) * xnew(O3 ) * xnew(NO )  &
        + rct(9,k) * xnew(NO ) * xnew(NO3 )  &
        + rct(9,k) * xnew(NO ) * xnew(NO3 )  &
        + rct(10,k) * xnew(NO ) * xnew(HO2 )  &
        + rct(11,k) * xnew(NO2 ) * xnew(NO3 )  &
        + rct(18,k) * xnew(OH ) * xnew(HONO )  &
        + rct(21,k) * xnew(CH3O2 ) * xnew(NO )  &
        + rct(30,k) * xnew(C2H5O2 ) * xnew(NO )  &
        + rct(34,k) * xnew(CH3COO2 ) * xnew(NO )  &
        + rct(40,k) * xnew(NO ) * xnew(SECC4H9O2 )  &
        + rct(40,k) * xnew(MEKO2 ) * xnew(NO )  &
        + rct(40,k) * xnew(ETRO2 ) * xnew(NO )  &
        + rct(40,k) * xnew(NO ) * xnew(PRRO2 )  &
        + rct(40,k) * xnew(OXYO2 ) * xnew(NO )  &
        + rct(40,k) * xnew(MALO2 ) * xnew(NO )  &
        + 0.86*rct(40,k) * xnew(ISRO2 ) * xnew(NO )  &
        + 0.95*rct(40,k) * xnew(MVKO2 ) * xnew(NO )  &
        + 0.95*rct(55,k) * xnew(MACRO2 ) * xnew(NO )  &
        + 2.5e-12 * xnew(MACRO2 ) * xnew(NO3 )  &
        + rct(57,k) * xnew(MACO3 ) * xnew(NO )  &
        + 1.9*rct(40,k) * xnew(ISNIR ) * xnew(NO )  &
        + 1.1*rct(40,k) * xnew(ISONO3 ) * xnew(NO )  &
        + 2.9e-11 * xnew(MPAN ) * xnew(OH )  &
        + rct(63,k) * xnew(OP ) * xnew(NO )  &
        + rct(65,k) * xnew(N2O5 )  &
        + rct(68,k) * xnew(PAN )  &
        + rct(68,k) * xnew(MPAN )  &
        + rcphot(IDHNO3,K) * xnew(HNO3 )  &
        + rcphot(IDNO3,K) * xnew(NO3 )  &
        + rcemis(NO2,k) 

      L =  &
          rct(6,k)* xnew(O3 )  &
        + rct(11,k)* xnew(NO3 )  &
        + rct(64,k)* xnew(NO3 )  &
        + rct(66,k)* xnew(OH )  &
        + rct(67,k)* xnew(CH3COO2 )  &
        + rct(67,k)* xnew(MACO3 )  &
        + rcphot(IDNO2,K) 

      xnew(NO2)=  ( xold(NO2) + dt2 * P) /(1.0 + dt2*L )

!-> PAN 

      P =   &
          rct(67,k) * xnew(CH3COO2 ) * xnew(NO2 ) 

      L =  &
          rct(68,k) 

      xnew(PAN)=  ( xold(PAN) + dt2 * P) /(1.0 + dt2*L )

!-> MPAN 

      P =   &
          rct(67,k) * xnew(MACO3 ) * xnew(NO2 ) 

      L =  &
          2.9e-11* xnew(OH )  &
        + rct(68,k) 

      xnew(MPAN)=  ( xold(MPAN) + dt2 * P) /(1.0 + dt2*L )

!-> NO3 

      P =   &
          rct(6,k) * xnew(O3 ) * xnew(NO2 )  &
        + rct(15,k) * xnew(OH ) * xnew(HNO3 )  &
        + rct(65,k) * xnew(N2O5 ) 

      L =  &
          rct(9,k)* xnew(NO )  &
        + rct(11,k)* xnew(NO2 )  &
        + rct(28,k)* xnew(HCHO )  &
        + 2.5e-12* xnew(MACRO2 )  &
        + rct(59,k)* xnew(C5H8 )  &
        + rct(64,k)* xnew(NO2 )  &
        + rcphot(IDNO3,K) 

      xnew(NO3)=  ( xold(NO3) + dt2 * P) /(1.0 + dt2*L )

!-> N2O5 

      P =   &
          rct(64,k) * xnew(NO2 ) * xnew(NO3 ) 

      L =  &
          rct(61,k)  &
        + rct(65,k) 

      xnew(N2O5)=  ( xold(N2O5) + dt2 * P) /(1.0 + dt2*L )

!-> ISONO3 

      P =   &
          rct(59,k) * xnew(C5H8 ) * xnew(NO3 )  &
        + 2.0e-11 * xnew(ISONO3H ) * xnew(OH ) 

      L =  &
          rct(40,k)* xnew(NO )  &
        + rct(47,k)* xnew(HO2 ) 

      xnew(ISONO3)=  ( xold(ISONO3) + dt2 * P) /(1.0 + dt2*L )

!-> HNO3 

      P =   &
          rct(28,k) * xnew(NO3 ) * xnew(HCHO )  &
        + 2.*rct(61,k) * xnew(N2O5 )  &
        + rct(66,k) * xnew(NO2 ) * xnew(OH ) 

      L =  &
          rct(15,k)* xnew(OH )  &
        + rct(62,k)  &
        + rcphot(IDHNO3,K) 

      xnew(HNO3)=  ( xold(HNO3) + dt2 * P) /(1.0 + dt2*L )

!-> HONO 

      P =   &
          rct(71,k) * xnew(OH ) * xnew(NO ) 

      L =  &
          rct(18,k)* xnew(OH )  &
        + 0.222*rcphot(IDNO2,K) 

      xnew(HONO)=  ( xold(HONO) + dt2 * P) /(1.0 + dt2*L )

!-> CH3COO2 

      P =   &
          0.95*rct(33,k) * xnew(OH ) * xnew(CH3CHO )  &
        + rct(32,k) * xnew(CH3COO2H ) * xnew(OH )  &
        + rct(40,k) * xnew(MEKO2 ) * xnew(NO )  &
        + rct(49,k) * xnew(OH ) * xnew(MGLYOX )  &
        + 0.684*rct(40,k) * xnew(MVKO2 ) * xnew(NO )  &
        + 0.41*rct(53,k) * xnew(MACR ) * xnew(O3 )  &
        + rct(57,k) * xnew(MACO3 ) * xnew(NO )  &
        + 1.51e-11 * xnew(MACO2H ) * xnew(OH )  &
        + rct(68,k) * xnew(PAN )  &
        + rcphot(IDRCOHCO,K) * xnew(MGLYOX )  &
        + rcphot(IDCH3O2H,K) * xnew(MEKO2H )  &
        + rcphot(IDCH3COX,K) * xnew(MEK ) 

      L =  &
          rct(34,k)* xnew(NO )  &
        + rct(35,k)* xnew(CH3O2 )  &
        + rct(36,k)* xnew(CH3COO2 )  &
        + rct(36,k)* xnew(CH3COO2 )  &
        + rct(37,k)* xnew(HO2 )  &
        + rct(67,k)* xnew(NO2 ) 

      xnew(CH3COO2)=  ( xold(CH3COO2) + dt2 * P) /(1.0 + dt2*L )

!-> MACR 

      P =   &
          0.67*rct(50,k) * xnew(C5H8 ) * xnew(O3 )  &
        + 0.32*rct(40,k) * xnew(ISRO2 ) * xnew(NO )  &
        + 0.1*rct(40,k) * xnew(ISONO3 ) * xnew(NO ) 

      L =  &
          rct(53,k)* xnew(O3 )  &
        + rct(54,k)* xnew(OH ) 

      xnew(MACR)=  ( xold(MACR) + dt2 * P) /(1.0 + dt2*L )

!-> ISNI 

      P =   &
          0.14*rct(40,k) * xnew(ISRO2 ) * xnew(NO )  &
        + 0.05*rct(40,k) * xnew(MVKO2 ) * xnew(NO )  &
        + 0.05*rct(55,k) * xnew(MACRO2 ) * xnew(NO )  &
        + 0.05*rct(40,k) * xnew(ISNIR ) * xnew(NO )  &
        + 0.85*rct(40,k) * xnew(ISONO3 ) * xnew(NO ) 

      L =  &
          5.96e-11* xnew(OH ) 

      xnew(ISNI)=  ( xold(ISNI) + dt2 * P) /(1.0 + dt2*L )

!-> ISNIR 

      P =   &
          5.96e-11 * xnew(ISNI ) * xnew(OH )  &
        + 3.7e-11 * xnew(ISNIRH ) * xnew(OH ) 

      L =  &
          rct(40,k)* xnew(NO )  &
        + rct(47,k)* xnew(HO2 ) 

      xnew(ISNIR)=  ( xold(ISNIR) + dt2 * P) /(1.0 + dt2*L )

!-> GLYOX 

      P =   &
          rct(40,k) * xnew(MALO2 ) * xnew(NO )  &
        + rcphot(IDCH3O2H,K) * xnew(MALO2H )  &
        + rcemis(GLYOX,k) 

      L =  &
          rct(48,k)* xnew(OH )  &
        + rcphot(IDHCOHCO,K) 

      xnew(GLYOX)=  ( xold(GLYOX) + dt2 * P) /(1.0 + dt2*L )

!-> MGLYOX 

      P =   &
          rct(40,k) * xnew(OXYO2 ) * xnew(NO )  &
        + rct(40,k) * xnew(MALO2 ) * xnew(NO )  &
        + 0.266*rct(40,k) * xnew(MVKO2 ) * xnew(NO )  &
        + 0.59*rct(53,k) * xnew(MACR ) * xnew(O3 )  &
        + rct(56,k) * xnew(ACETOL ) * xnew(OH )  &
        + 0.82*rct(58,k) * xnew(MVK ) * xnew(O3 )  &
        + rcphot(IDCH3O2H,K) * xnew(OXYO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(MALO2H )  &
        + rcemis(MGLYOX,k) 

      L =  &
          rct(49,k)* xnew(OH )  &
        + rcphot(IDRCOHCO,K) 

      xnew(MGLYOX)=  ( xold(MGLYOX) + dt2 * P) /(1.0 + dt2*L )

!-> MAL 

      P =   &
          rct(40,k) * xnew(OXYO2 ) * xnew(NO )  &
        + rcphot(IDCH3O2H,K) * xnew(OXYO2H ) 

      L =  &
          5.58e-11* xnew(OH ) 

      xnew(MAL)=  ( xold(MAL) + dt2 * P) /(1.0 + dt2*L )

!-> MEK 

      P =   &
          0.65*rct(40,k) * xnew(NO ) * xnew(SECC4H9O2 )  &
        + 2.15e-11 * xnew(BURO2H ) * xnew(OH )  &
        + 0.65*rcphot(IDCH3O2H,K) * xnew(BURO2H )  &
        + rcemis(MEK,k) 

      L =  &
          rct(42,k)* xnew(OH )  &
        + rcphot(IDCH3COX,K) 

      xnew(MEK)=  ( xold(MEK) + dt2 * P) /(1.0 + dt2*L )

!-> MVK 

      P =   &
          0.26*rct(50,k) * xnew(C5H8 ) * xnew(O3 )  &
        + 0.42*rct(40,k) * xnew(ISRO2 ) * xnew(NO )  &
        + 0.05*rct(40,k) * xnew(ISONO3 ) * xnew(NO ) 

      L =  &
          rct(52,k)* xnew(OH )  &
        + rct(58,k)* xnew(O3 ) 

      xnew(MVK)=  ( xold(MVK) + dt2 * P) /(1.0 + dt2*L )

!-> HCHO 

      P =   &
          rct(21,k) * xnew(CH3O2 ) * xnew(NO )  &
        + 2.*rct(22,k) * xnew(CH3O2 ) * xnew(CH3O2 )  &
        + rct(23,k) * xnew(CH3O2 ) * xnew(CH3O2 )  &
        + rct(24,k) * xnew(OH ) * xnew(CH3OH )  &
        + 0.1*rct(25,k) * xnew(HO2 ) * xnew(CH3O2 )  &
        + 0.4*rct(26,k) * xnew(CH3O2H ) * xnew(OH )  &
        + rct(35,k) * xnew(CH3O2 ) * xnew(CH3COO2 )  &
        + 2.*rct(40,k) * xnew(ETRO2 ) * xnew(NO )  &
        + 1.14*rct(43,k) * xnew(C2H4 ) * xnew(O3 )  &
        + 0.545*rct(44,k) * xnew(O3 ) * xnew(C3H6 )  &
        + rct(40,k) * xnew(NO ) * xnew(PRRO2 )  &
        + 0.8*rct(50,k) * xnew(C5H8 ) * xnew(O3 )  &
        + 0.74*rct(40,k) * xnew(ISRO2 ) * xnew(NO )  &
        + 0.266*rct(40,k) * xnew(MVKO2 ) * xnew(NO )  &
        + 0.534*rct(53,k) * xnew(MACR ) * xnew(O3 )  &
        + 2.5e-12 * xnew(MACRO2 ) * xnew(NO3 )  &
        + rct(57,k) * xnew(MACO3 ) * xnew(NO )  &
        + 0.8*rct(58,k) * xnew(MVK ) * xnew(O3 )  &
        + 0.15*rct(40,k) * xnew(ISONO3 ) * xnew(NO )  &
        + 1.51e-11 * xnew(MACO2H ) * xnew(OH )  &
        + rcphot(IDCH3O2H,K) * xnew(CH3O2H )  &
        + 0.1*rcphot(IDHCOHCO,K) * xnew(GLYOX )  &
        + 1.56*rcphot(IDCH3O2H,K) * xnew(ETRO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(PRRO2H )  &
        + rcemis(HCHO,k) 

      L =  &
          rct(27,k)* xnew(OH )  &
        + rct(28,k)* xnew(NO3 )  &
        + rcphot(IDACH2O,K)  &
        + rcphot(IDBCH2O,K) 

      xnew(HCHO)=  ( xold(HCHO) + dt2 * P) /(1.0 + dt2*L )

!-> CH3CHO 

      P =   &
          rct(30,k) * xnew(C2H5O2 ) * xnew(NO )  &
        + 8.01e-12 * xnew(C2H5OOH ) * xnew(OH )  &
        + rct(38,k) * xnew(OH ) * xnew(C2H5OH )  &
        + 0.35*rct(40,k) * xnew(NO ) * xnew(SECC4H9O2 )  &
        + rct(40,k) * xnew(MEKO2 ) * xnew(NO )  &
        + 1.38e-11 * xnew(ETRO2H ) * xnew(OH )  &
        + 0.545*rct(44,k) * xnew(O3 ) * xnew(C3H6 )  &
        + rct(40,k) * xnew(NO ) * xnew(PRRO2 )  &
        + 0.684*rct(40,k) * xnew(MVKO2 ) * xnew(NO )  &
        + 0.04*rct(58,k) * xnew(MVK ) * xnew(O3 )  &
        + 0.95*rct(40,k) * xnew(ISNIR ) * xnew(NO )  &
        + rcphot(IDCH3O2H,K) * xnew(C2H5OOH )  &
        + 0.22*rcphot(IDCH3O2H,K) * xnew(ETRO2H )  &
        + 0.35*rcphot(IDCH3O2H,K) * xnew(BURO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(PRRO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(MEKO2H )  &
        + rcemis(CH3CHO,k) 

      L =  &
          rct(33,k)* xnew(OH )  &
        + rcphot(IDCH3CHO,K) 

      xnew(CH3CHO)=  ( xold(CH3CHO) + dt2 * P) /(1.0 + dt2*L )

!-> C2H6 

      P =   &
          rcemis(C2H6,k) 

      L =  &
          rct(29,k)* xnew(OH ) 

      xnew(C2H6)=  ( xold(C2H6) + dt2 * P) /(1.0 + dt2*L )

!-> NC4H10 

      P =   &
          rcemis(NC4H10,k) 

      L =  &
          rct(39,k)* xnew(OH ) 

      xnew(NC4H10)=  ( xold(NC4H10) + dt2 * P) /(1.0 + dt2*L )

!-> C2H4 

      P =   &
          rcemis(C2H4,k) 

      L =  &
          rct(43,k)* xnew(O3 )  &
        + rct(69,k)* xnew(OH ) 

      xnew(C2H4)=  ( xold(C2H4) + dt2 * P) /(1.0 + dt2*L )

!-> C3H6 

      P =   &
          0.07*rct(50,k) * xnew(C5H8 ) * xnew(O3 )  &
        + rcemis(C3H6,k) 

      L =  &
          rct(44,k)* xnew(O3 )  &
        + rct(70,k)* xnew(OH ) 

      xnew(C3H6)=  ( xold(C3H6) + dt2 * P) /(1.0 + dt2*L )

!-> OXYL 

      P =   &
          rcemis(OXYL,k) 

      L =  &
          1.36e-11* xnew(OH ) 

      xnew(OXYL)=  ( xold(OXYL) + dt2 * P) /(1.0 + dt2*L )

!-> C5H8 

      P =   &
          rcemis(C5H8,k)  &
        + 0 !Skip bio rate since rcemis exists 

      L =  &
          rct(50,k)* xnew(O3 )  &
        + rct(51,k)* xnew(OH )  &
        + rct(59,k)* xnew(NO3 ) 

      xnew(C5H8)=  ( xold(C5H8) + dt2 * P) /(1.0 + dt2*L )

!-> APINENE 

      P =   &
          rcemis(APINENE,k) 

      L =  &
          rct(72,k)*xnew(O3)  &
        + rct(73,k)*xnew(OH)  &
        + rct(74,k)*xnew(NO3) 

      xnew(APINENE)=  ( xold(APINENE) + dt2 * P) /(1.0 + dt2*L )

!-> CH3O2H 

      P =   &
          0.9*rct(25,k) * xnew(HO2 ) * xnew(CH3O2 ) 

      L =  &
          rct(26,k)* xnew(OH )  &
        + 1.0e-5  &
        + rcphot(IDCH3O2H,K) 

      xnew(CH3O2H)=  ( xold(CH3O2H) + dt2 * P) /(1.0 + dt2*L )

!-> C2H5OOH 

      P =   &
          rct(31,k) * xnew(C2H5O2 ) * xnew(HO2 ) 

      L =  &
          8.01e-12* xnew(OH )  &
        + rct(32,k)* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(C2H5OOH)=  ( xold(C2H5OOH) + dt2 * P) /(1.0 + dt2*L )

!-> BURO2H 

      P =   &
          0.95*rct(41,k) * xnew(SECC4H9O2 ) * xnew(HO2 ) 

      L =  &
          rct(32,k)* xnew(OH )  &
        + 2.15e-11* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(BURO2H)=  ( xold(BURO2H) + dt2 * P) /(1.0 + dt2*L )

!-> ETRO2H 

      P =   &
          1.2e-11 * xnew(ETRO2 ) * xnew(HO2 ) 

      L =  &
          1.38e-11* xnew(OH )  &
        + rct(32,k)* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(ETRO2H)=  ( xold(ETRO2H) + dt2 * P) /(1.0 + dt2*L )

!-> PRRO2H 

      P =   &
          0.795*rct(45,k) * xnew(PRRO2 ) * xnew(HO2 ) 

      L =  &
          2.44e-11* xnew(OH )  &
        + rct(32,k)* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(PRRO2H)=  ( xold(PRRO2H) + dt2 * P) /(1.0 + dt2*L )

!-> OXYO2H 

      P =   &
          0.33*rct(46,k) * xnew(OXYO2 ) * xnew(HO2 ) 

      L =  &
          4.2e-11* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(OXYO2H)=  ( xold(OXYO2H) + dt2 * P) /(1.0 + dt2*L )

!-> MEKO2H 

      P =   &
          rct(41,k) * xnew(MEKO2 ) * xnew(HO2 ) 

      L =  &
          rct(32,k)* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(MEKO2H)=  ( xold(MEKO2H) + dt2 * P) /(1.0 + dt2*L )

!-> MALO2H 

      P =   &
          rct(47,k) * xnew(MALO2 ) * xnew(HO2 ) 

      L =  &
          rct(32,k)* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(MALO2H)=  ( xold(MALO2H) + dt2 * P) /(1.0 + dt2*L )

!-> MVKO2H 

      P =   &
          rct(41,k) * xnew(MVKO2 ) * xnew(HO2 ) 

      L =  &
          2.2e-11* xnew(OH ) 

      xnew(MVKO2H)=  ( xold(MVKO2H) + dt2 * P) /(1.0 + dt2*L )

!-> MACROOH 

      P =   &
          rct(41,k) * xnew(MACRO2 ) * xnew(HO2 ) 

      L =  &
          2.82e-11* xnew(OH ) 

      xnew(MACROOH)=  ( xold(MACROOH) + dt2 * P) /(1.0 + dt2*L )

!-> MACO3H 

      P =   &
          0.71*rct(60,k) * xnew(MACO3 ) * xnew(HO2 ) 

      L =  &
          1.87e-11* xnew(OH ) 

      xnew(MACO3H)=  ( xold(MACO3H) + dt2 * P) /(1.0 + dt2*L )

!-> MACO2H 

      P =   &
          0.29*rct(60,k) * xnew(MACO3 ) * xnew(HO2 ) 

      L =  &
          1.51e-11* xnew(OH ) 

      xnew(MACO2H)=  ( xold(MACO2H) + dt2 * P) /(1.0 + dt2*L )

!-> ISRO2H 

      P =   &
          0.97*rct(47,k) * xnew(ISRO2 ) * xnew(HO2 ) 

      L =  &
          7.5e-11* xnew(OH ) 

      xnew(ISRO2H)=  ( xold(ISRO2H) + dt2 * P) /(1.0 + dt2*L )

!-> H2O2 

      P =   &
          rct(16,k) * xnew(HO2 ) * xnew(HO2 )  &
        + rct(17,k) * xnew(HO2 ) * xnew(HO2 )  &
        + 0.14*rct(43,k) * xnew(C2H4 ) * xnew(O3 )  &
        + 0.09*rct(44,k) * xnew(O3 ) * xnew(C3H6 )  &
        + 0.124*rct(53,k) * xnew(MACR ) * xnew(O3 ) 

      L =  &
          rct(13,k)* xnew(OH )  &
        + AQRCK(ICLRC1,K)* xnew(SO2 )  &
        + 1.0e-5  &
        + rcphot(IDH2O2,K) 

      xnew(H2O2)=  ( xold(H2O2) + dt2 * P) /(1.0 + dt2*L )

!-> CH3COO2H 

      P =   &
          0.41*rct(37,k) * xnew(CH3COO2 ) * xnew(HO2 ) 

      L =  &
          rct(32,k)* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(CH3COO2H)=  ( xold(CH3COO2H) + dt2 * P) /(1.0 + dt2*L )

!-> ISONO3H 

      P =   &
          rct(47,k) * xnew(ISONO3 ) * xnew(HO2 ) 

      L =  &
          2.0e-11* xnew(OH ) 

      xnew(ISONO3H)=  ( xold(ISONO3H) + dt2 * P) /(1.0 + dt2*L )

!-> ISNIRH 

      P =   &
          rct(47,k) * xnew(ISNIR ) * xnew(HO2 ) 

      L =  &
          3.7e-11* xnew(OH ) 

      xnew(ISNIRH)=  ( xold(ISNIRH) + dt2 * P) /(1.0 + dt2*L )

!-> CH3OH 

      P =   &
          rct(23,k) * xnew(CH3O2 ) * xnew(CH3O2 )  &
        + rcemis(CH3OH,k) 

      L =  &
          rct(24,k)* xnew(OH ) 

      xnew(CH3OH)=  ( xold(CH3OH) + dt2 * P) /(1.0 + dt2*L )

!-> C2H5OH 

      P =   &
          rcemis(C2H5OH,k) 

      L =  &
          rct(38,k)* xnew(OH ) 

      xnew(C2H5OH)=  ( xold(C2H5OH) + dt2 * P) /(1.0 + dt2*L )

!-> ACETOL 

      P =   &
          2.44e-11 * xnew(PRRO2H ) * xnew(OH )  &
        + 0.95*rct(55,k) * xnew(MACRO2 ) * xnew(NO )  &
        + 2.5e-12 * xnew(MACRO2 ) * xnew(NO3 )  &
        + 0.95*rct(40,k) * xnew(ISNIR ) * xnew(NO )  &
        + 2.9e-11 * xnew(MPAN ) * xnew(OH ) 

      L =  &
          rct(56,k)* xnew(OH ) 

      xnew(ACETOL)=  ( xold(ACETOL) + dt2 * P) /(1.0 + dt2*L )

!-> H2 

      P =   &
          rcphot(IDBCH2O,K) * xnew(HCHO ) 

      L =  &
          rct(14,k)* xnew(OH ) 

      xnew(H2)=  ( xold(H2) + dt2 * P) /(1.0 + dt2*L )

!-> CO 

      P =   &
          rct(27,k) * xnew(OH ) * xnew(HCHO )  &
        + rct(28,k) * xnew(NO3 ) * xnew(HCHO )  &
        + 0.05*rct(33,k) * xnew(OH ) * xnew(CH3CHO )  &
        + 0.63*rct(43,k) * xnew(C2H4 ) * xnew(O3 )  &
        + 0.56*rct(44,k) * xnew(O3 ) * xnew(C3H6 )  &
        + rct(48,k) * xnew(OH ) * xnew(GLYOX )  &
        + rct(48,k) * xnew(OH ) * xnew(GLYOX )  &
        + rct(49,k) * xnew(OH ) * xnew(MGLYOX )  &
        + 0.05*rct(50,k) * xnew(C5H8 ) * xnew(O3 )  &
        + 0.82*rct(53,k) * xnew(MACR ) * xnew(O3 )  &
        + 0.95*rct(55,k) * xnew(MACRO2 ) * xnew(NO )  &
        + 0.05*rct(58,k) * xnew(MVK ) * xnew(O3 )  &
        + 2.9e-11 * xnew(MPAN ) * xnew(OH )  &
        + rcphot(IDACH2O,K) * xnew(HCHO )  &
        + rcphot(IDBCH2O,K) * xnew(HCHO )  &
        + rcphot(IDCH3CHO,K) * xnew(CH3CHO )  &
        + 1.9*rcphot(IDHCOHCO,K) * xnew(GLYOX )  &
        + rcphot(IDRCOHCO,K) * xnew(MGLYOX )  &
        + rcemis(CO,k) 

      L =  &
          rct(20,k)* xnew(OH ) 

      xnew(CO)=  ( xold(CO) + dt2 * P) /(1.0 + dt2*L )

!-> CH4 

      P =   &
          0.1*rct(44,k) * xnew(O3 ) * xnew(C3H6 ) 

      L =  &
          rct(19,k)* xnew(OH ) 

      xnew(CH4)=  ( xold(CH4) + dt2 * P) /(1.0 + dt2*L )

!-> SO2 

      P =   &
          rcemis(SO2,k) 

      L =  &
          2e-12*AQRCK(ICLOHSO2,K)* xnew(OH )  &
        + AQRCK(ICLRC1,K)* xnew(H2O2 )  &
        + AQRCK(ICLRC2,K)* xnew(O3 )  &
        + AQRCK(ICLRC3,K) 

      xnew(SO2)=  ( xold(SO2) + dt2 * P) /(1.0 + dt2*L )
