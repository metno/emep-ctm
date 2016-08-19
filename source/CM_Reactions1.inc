
!-> OD 

      P =   &
          rcphot(IDBO3,K) * xnew(O3 ) 

      L =  &
          rcmisc(2,k)  &
        + rcmisc(3,k)  &
        + rcmisc(4,k) 

      xnew(OD)=  max(0.0, ( xold(OD) + dt2 * P)) /(1.0 + dt2*L )

!-> OP 

      P =   &
          rcmisc(2,k) * xnew(OD )  &
        + rcmisc(3,k) * xnew(OD )  &
        + 0.3*rct(34,k) * xnew(C5H8 ) * xnew(O3 )  &
        + 0.2*rct(42,k) * xnew(MVK ) * xnew(O3 )  &
        + rcphot(IDAO3,K) * xnew(O3 )  &
        + rcphot(IDNO2,K) * xnew(NO2 )  &
        + rcphot(IDNO3,K) * xnew(NO3 ) 

      L =  &
          rcmisc(1,k)  &
        + rcmisc(20,k)* xnew(NO ) 

      xnew(OP)=  max(0.0, ( xold(OP) + dt2 * P)) /(1.0 + dt2*L )

!-> OH 

      P =   &
          2.*rcmisc(4,k) * xnew(OD )  &
        + rcmisc(5,k) * xnew(O3 ) * xnew(HO2 )  &
        + rct(5,k) * xnew(NO ) * xnew(HO2 )  &
        + 0.4*rct(15,k) * xnew(CH3O2H ) * xnew(OH )  &
        + 8.01e-12 * xnew(C2H5OOH ) * xnew(OH )  &
        + 0.44*rct(25,k) * xnew(CH3COO2 ) * xnew(HO2 )  &
        + 2.15e-11 * xnew(BURO2H ) * xnew(OH )  &
        + 1.38e-11 * xnew(ETRO2H ) * xnew(OH )  &
        + 0.13*rct(28,k) * xnew(C2H4 ) * xnew(O3 )  &
        + 0.36*rct(29,k) * xnew(O3 ) * xnew(C3H6 )  &
        + 2.44e-11 * xnew(PRRO2H ) * xnew(OH )  &
        + 0.55*rct(34,k) * xnew(C5H8 ) * xnew(O3 )  &
        + 7.5e-11 * xnew(ISRO2H ) * xnew(OH )  &
        + 0.82*rct(37,k) * xnew(MACR ) * xnew(O3 )  &
        + 0.08*rct(42,k) * xnew(MVK ) * xnew(O3 )  &
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
          rct(3,k)* xnew(O3 )  &
        + rct(7,k)* xnew(HO2 )  &
        + rct(8,k)* xnew(H2O2 )  &
        + rct(9,k)* xnew(H2 )  &
        + rcmisc(6,k)* xnew(HNO3 )  &
        + rct(10,k)* xnew(HONO )  &
        + 2e-12*AQRCK(ICLOHSO2,K)* xnew(SO2 )  &
        + rcmisc(9,k)* xnew(CH4 )  &
        + rcmisc(10,k)* xnew(CO )  &
        + rcmisc(11,k)* xnew(CH3OH )  &
        + rct(15,k)* xnew(CH3O2H )  &
        + rcmisc(12,k)* xnew(HCHO )  &
        + rct(17,k)* xnew(C2H6 )  &
        + 8.01e-12* xnew(C2H5OOH )  &
        + rct(20,k)* xnew(C2H5OOH )  &
        + rct(21,k)* xnew(CH3CHO )  &
        + rct(20,k)* xnew(CH3COO2H )  &
        + rcmisc(13,k)* xnew(C2H5OH )  &
        + rcmisc(14,k)* xnew(NC4H10 )  &
        + rcmisc(15,k)* xnew(MEK )  &
        + rct(20,k)* xnew(MEKO2H )  &
        + rct(20,k)* xnew(BURO2H )  &
        + 2.15e-11* xnew(BURO2H )  &
        + 1.38e-11* xnew(ETRO2H )  &
        + rct(20,k)* xnew(ETRO2H )  &
        + 2.44e-11* xnew(PRRO2H )  &
        + rct(20,k)* xnew(PRRO2H )  &
        + 1.36e-11* xnew(OXYL )  &
        + 4.2e-11* xnew(OXYO2H )  &
        + 5.58e-11* xnew(MAL )  &
        + rct(20,k)* xnew(MALO2H )  &
        + rcmisc(16,k)* xnew(GLYOX )  &
        + rct(33,k)* xnew(MGLYOX )  &
        + rct(35,k)* xnew(C5H8 )  &
        + rct(36,k)* xnew(MVK )  &
        + 7.5e-11* xnew(ISRO2H )  &
        + rct(38,k)* xnew(MACR )  &
        + 2.82e-11* xnew(MACROOH )  &
        + rct(40,k)* xnew(ACETOL )  &
        + 5.96e-11* xnew(ISNI )  &
        + 1.87e-11* xnew(MACO3H )  &
        + 1.51e-11* xnew(MACO2H )  &
        + 2.0e-11* xnew(ISONO3H )  &
        + 2.2e-11* xnew(MVKO2H )  &
        + 3.7e-11* xnew(ISNIRH )  &
        + 2.9e-11* xnew(MPAN )  &
        + rcmisc(23,k)* xnew(NO2 )  &
        + rcmisc(26,k)* xnew(C2H4 )  &
        + rcmisc(27,k)* xnew(C3H6 )  &
        + rcmisc(30,k)* xnew(NO ) 

      xnew(OH)=  max(0.0, ( xold(OH) + dt2 * P)) /(1.0 + dt2*L )

!-> HO2 

      P =   &
          rct(3,k) * xnew(O3 ) * xnew(OH )  &
        + rct(8,k) * xnew(OH ) * xnew(H2O2 )  &
        + rct(9,k) * xnew(OH ) * xnew(H2 )  &
        + 2e-12*AQRCK(ICLOHSO2,K) * xnew(OH ) * xnew(SO2 )  &
        + rcmisc(10,k) * xnew(OH ) * xnew(CO )  &
        + rct(11,k) * xnew(CH3O2 ) * xnew(NO )  &
        + 2.*rct(12,k) * xnew(CH3O2 ) * xnew(CH3O2 )  &
        + rcmisc(11,k) * xnew(OH ) * xnew(CH3OH )  &
        + rcmisc(12,k) * xnew(OH ) * xnew(HCHO )  &
        + rct(16,k) * xnew(NO3 ) * xnew(HCHO )  &
        + rct(18,k) * xnew(C2H5O2 ) * xnew(NO )  &
        + 0.9*rct(23,k) * xnew(CH3O2 ) * xnew(CH3COO2 )  &
        + rcmisc(13,k) * xnew(OH ) * xnew(C2H5OH )  &
        + 0.65*rct(26,k) * xnew(NO ) * xnew(SECC4H9O2 )  &
        + rct(26,k) * xnew(ETRO2 ) * xnew(NO )  &
        + 0.13*rct(28,k) * xnew(C2H4 ) * xnew(O3 )  &
        + 0.28*rct(29,k) * xnew(O3 ) * xnew(C3H6 )  &
        + rct(26,k) * xnew(NO ) * xnew(PRRO2 )  &
        + rct(26,k) * xnew(OXYO2 ) * xnew(NO )  &
        + rct(26,k) * xnew(MALO2 ) * xnew(NO )  &
        + rcmisc(16,k) * xnew(OH ) * xnew(GLYOX )  &
        + 0.06*rct(34,k) * xnew(C5H8 ) * xnew(O3 )  &
        + 0.78*rct(26,k) * xnew(ISRO2 ) * xnew(NO )  &
        + 0.95*rct(26,k) * xnew(MVKO2 ) * xnew(NO )  &
        + 0.41*rct(37,k) * xnew(MACR ) * xnew(O3 )  &
        + 0.95*rct(39,k) * xnew(MACRO2 ) * xnew(NO )  &
        + 2.5e-12 * xnew(MACRO2 ) * xnew(NO3 )  &
        + rct(40,k) * xnew(ACETOL ) * xnew(OH )  &
        + 0.06*rct(42,k) * xnew(MVK ) * xnew(O3 )  &
        + 0.05*rct(26,k) * xnew(ISNIR ) * xnew(NO )  &
        + 0.8*rct(26,k) * xnew(ISONO3 ) * xnew(NO )  &
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
          rcmisc(5,k)* xnew(O3 )  &
        + rct(5,k)* xnew(NO )  &
        + rct(7,k)* xnew(OH )  &
        + rcmisc(7,k)* xnew(HO2 )  &
        + rcmisc(7,k)* xnew(HO2 )  &
        + rcmisc(8,k)* xnew(HO2 )  &
        + rcmisc(8,k)* xnew(HO2 )  &
        + rct(14,k)* xnew(CH3O2 )  &
        + rct(19,k)* xnew(C2H5O2 )  &
        + rct(25,k)* xnew(CH3COO2 )  &
        + rct(27,k)* xnew(SECC4H9O2 )  &
        + rct(27,k)* xnew(MEKO2 )  &
        + 1.2e-11* xnew(ETRO2 )  &
        + rct(30,k)* xnew(PRRO2 )  &
        + rct(31,k)* xnew(OXYO2 )  &
        + rct(32,k)* xnew(MALO2 )  &
        + rct(32,k)* xnew(ISRO2 )  &
        + rct(27,k)* xnew(MACRO2 )  &
        + rct(32,k)* xnew(ISONO3 )  &
        + rct(27,k)* xnew(MVKO2 )  &
        + rct(44,k)* xnew(MACO3 )  &
        + rct(32,k)* xnew(ISNIR ) 

      xnew(HO2)=  max(0.0, ( xold(HO2) + dt2 * P)) /(1.0 + dt2*L )

!-> CH3O2 

      P =   &
          rcmisc(9,k) * xnew(OH ) * xnew(CH4 )  &
        + 0.6*rct(15,k) * xnew(CH3O2H ) * xnew(OH )  &
        + 0.05*rct(21,k) * xnew(OH ) * xnew(CH3CHO )  &
        + rct(22,k) * xnew(CH3COO2 ) * xnew(NO )  &
        + 0.9*rct(23,k) * xnew(CH3O2 ) * xnew(CH3COO2 )  &
        + rct(24,k) * xnew(CH3COO2 ) * xnew(CH3COO2 )  &
        + rct(24,k) * xnew(CH3COO2 ) * xnew(CH3COO2 )  &
        + 0.44*rct(25,k) * xnew(CH3COO2 ) * xnew(HO2 )  &
        + 0.28*rct(29,k) * xnew(O3 ) * xnew(C3H6 )  &
        + rcphot(IDCH3CHO,K) * xnew(CH3CHO )  &
        + rcphot(IDCH3O2H,K) * xnew(CH3COO2H ) 

      L =  &
          rct(11,k)* xnew(NO )  &
        + rct(12,k)* xnew(CH3O2 )  &
        + rct(12,k)* xnew(CH3O2 )  &
        + rct(13,k)* xnew(CH3O2 )  &
        + rct(13,k)* xnew(CH3O2 )  &
        + rct(14,k)* xnew(HO2 )  &
        + rct(23,k)* xnew(CH3COO2 ) 

      xnew(CH3O2)=  max(0.0, ( xold(CH3O2) + dt2 * P)) /(1.0 + dt2*L )

!-> C2H5O2 

      P =   &
          rct(17,k) * xnew(OH ) * xnew(C2H6 )  &
        + rct(20,k) * xnew(C2H5OOH ) * xnew(OH )  &
        + 0.35*rct(26,k) * xnew(NO ) * xnew(SECC4H9O2 )  &
        + 0.35*rcphot(IDCH3O2H,K) * xnew(BURO2H )  &
        + rcphot(IDCH3COX,K) * xnew(MEK ) 

      L =  &
          rct(18,k)* xnew(NO )  &
        + rct(19,k)* xnew(HO2 ) 

      xnew(C2H5O2)=  max(0.0, ( xold(C2H5O2) + dt2 * P)) /(1.0 + dt2*L )

!-> SECC4H9O2 

      P =   &
          rcmisc(14,k) * xnew(OH ) * xnew(NC4H10 )  &
        + rct(20,k) * xnew(BURO2H ) * xnew(OH ) 

      L =  &
          rct(26,k)* xnew(NO )  &
        + rct(27,k)* xnew(HO2 ) 

      xnew(SECC4H9O2)=  max(0.0, ( xold(SECC4H9O2) + dt2 * P)) /(1.0 + dt2*L )

!-> ISRO2 

      P =   &
          rct(35,k) * xnew(C5H8 ) * xnew(OH )  &
        + 0.12*rct(26,k) * xnew(ISRO2 ) * xnew(NO )  &
        + 7.5e-11 * xnew(ISRO2H ) * xnew(OH ) 

      L =  &
          rct(26,k)* xnew(NO )  &
        + rct(32,k)* xnew(HO2 ) 

      xnew(ISRO2)=  max(0.0, ( xold(ISRO2) + dt2 * P)) /(1.0 + dt2*L )

!-> ETRO2 

      P =   &
          rct(20,k) * xnew(ETRO2H ) * xnew(OH )  &
        + rcmisc(26,k) * xnew(C2H4 ) * xnew(OH ) 

      L =  &
          rct(26,k)* xnew(NO )  &
        + 1.2e-11* xnew(HO2 ) 

      xnew(ETRO2)=  max(0.0, ( xold(ETRO2) + dt2 * P)) /(1.0 + dt2*L )

!-> PRRO2 

      P =   &
          rct(20,k) * xnew(PRRO2H ) * xnew(OH )  &
        + rcmisc(27,k) * xnew(OH ) * xnew(C3H6 ) 

      L =  &
          rct(26,k)* xnew(NO )  &
        + rct(30,k)* xnew(HO2 ) 

      xnew(PRRO2)=  max(0.0, ( xold(PRRO2) + dt2 * P)) /(1.0 + dt2*L )

!-> OXYO2 

      P =   &
          1.36e-11 * xnew(OXYL ) * xnew(OH )  &
        + 4.2e-11 * xnew(OXYO2H ) * xnew(OH ) 

      L =  &
          rct(26,k)* xnew(NO )  &
        + rct(31,k)* xnew(HO2 ) 

      xnew(OXYO2)=  max(0.0, ( xold(OXYO2) + dt2 * P)) /(1.0 + dt2*L )

!-> MEKO2 

      P =   &
          rcmisc(15,k) * xnew(OH ) * xnew(MEK )  &
        + rct(20,k) * xnew(MEKO2H ) * xnew(OH ) 

      L =  &
          rct(26,k)* xnew(NO )  &
        + rct(27,k)* xnew(HO2 ) 

      xnew(MEKO2)=  max(0.0, ( xold(MEKO2) + dt2 * P)) /(1.0 + dt2*L )

!-> MALO2 

      P =   &
          5.58e-11 * xnew(MAL ) * xnew(OH )  &
        + rct(20,k) * xnew(MALO2H ) * xnew(OH ) 

      L =  &
          rct(26,k)* xnew(NO )  &
        + rct(32,k)* xnew(HO2 ) 

      xnew(MALO2)=  max(0.0, ( xold(MALO2) + dt2 * P)) /(1.0 + dt2*L )

!-> MVKO2 

      P =   &
          rct(36,k) * xnew(MVK ) * xnew(OH )  &
        + 2.2e-11 * xnew(MVKO2H ) * xnew(OH ) 

      L =  &
          rct(26,k)* xnew(NO )  &
        + rct(27,k)* xnew(HO2 ) 

      xnew(MVKO2)=  max(0.0, ( xold(MVKO2) + dt2 * P)) /(1.0 + dt2*L )

!-> MACRO2 

      P =   &
          0.5*rct(38,k) * xnew(MACR ) * xnew(OH )  &
        + 2.82e-11 * xnew(MACROOH ) * xnew(OH ) 

      L =  &
          rct(39,k)* xnew(NO )  &
        + 2.5e-12* xnew(NO3 )  &
        + rct(27,k)* xnew(HO2 ) 

      xnew(MACRO2)=  max(0.0, ( xold(MACRO2) + dt2 * P)) /(1.0 + dt2*L )

!-> MACO3 

      P =   &
          0.5*rct(38,k) * xnew(MACR ) * xnew(OH )  &
        + 1.87e-11 * xnew(MACO3H ) * xnew(OH )  &
        + rcmisc(29,k) * xnew(MPAN ) 

      L =  &
          rct(41,k)* xnew(NO )  &
        + rct(44,k)* xnew(HO2 )  &
        + rcmisc(28,k)* xnew(NO2 ) 

      xnew(MACO3)=  max(0.0, ( xold(MACO3) + dt2 * P)) /(1.0 + dt2*L )

!-> O3 

      P =   &
          rcmisc(1,k) * xnew(OP )  &
        + 0.15*rct(25,k) * xnew(CH3COO2 ) * xnew(HO2 )  &
        + 0.29*rct(44,k) * xnew(MACO3 ) * xnew(HO2 ) 

      L =  &
          rct(1,k)* xnew(NO )  &
        + rct(2,k)* xnew(NO2 )  &
        + rct(3,k)* xnew(OH )  &
        + rcmisc(5,k)* xnew(HO2 )  &
        + AQRCK(ICLRC2,K)* xnew(SO2 )  &
        + rct(28,k)* xnew(C2H4 )  &
        + rct(29,k)* xnew(C3H6 )  &
        + rct(34,k)* xnew(C5H8 )  &
        + rct(37,k)* xnew(MACR )  &
        + rct(42,k)* xnew(MVK )  &
        + rcphot(IDAO3,K)  &
        + rcphot(IDBO3,K) 

      xnew(O3)=  max(0.0, ( xold(O3) + dt2 * P)) /(1.0 + dt2*L )

!-> NO 

      P =   &
          rct(6,k) * xnew(NO2 ) * xnew(NO3 )  &
        + rcphot(IDNO2,K) * xnew(NO2 )  &
        + 0.222*rcphot(IDNO2,K) * xnew(HONO )  &
        + rcemis(NO,k) 

      L =  &
          rct(1,k)* xnew(O3 )  &
        + rct(4,k)* xnew(NO3 )  &
        + rct(5,k)* xnew(HO2 )  &
        + rct(11,k)* xnew(CH3O2 )  &
        + rct(18,k)* xnew(C2H5O2 )  &
        + rct(22,k)* xnew(CH3COO2 )  &
        + rct(26,k)* xnew(SECC4H9O2 )  &
        + rct(26,k)* xnew(MEKO2 )  &
        + rct(26,k)* xnew(ETRO2 )  &
        + rct(26,k)* xnew(PRRO2 )  &
        + rct(26,k)* xnew(OXYO2 )  &
        + rct(26,k)* xnew(MALO2 )  &
        + rct(26,k)* xnew(ISRO2 )  &
        + rct(26,k)* xnew(MVKO2 )  &
        + rct(39,k)* xnew(MACRO2 )  &
        + rct(41,k)* xnew(MACO3 )  &
        + rct(26,k)* xnew(ISNIR )  &
        + rct(26,k)* xnew(ISONO3 )  &
        + rcmisc(20,k)* xnew(OP )  &
        + rcmisc(30,k)* xnew(OH ) 

      xnew(NO)=  max(0.0, ( xold(NO) + dt2 * P)) /(1.0 + dt2*L )

!-> NO2 

      P =   &
          rct(1,k) * xnew(O3 ) * xnew(NO )  &
        + rct(4,k) * xnew(NO ) * xnew(NO3 )  &
        + rct(4,k) * xnew(NO ) * xnew(NO3 )  &
        + rct(5,k) * xnew(NO ) * xnew(HO2 )  &
        + rct(6,k) * xnew(NO2 ) * xnew(NO3 )  &
        + rct(10,k) * xnew(OH ) * xnew(HONO )  &
        + rct(11,k) * xnew(CH3O2 ) * xnew(NO )  &
        + rct(18,k) * xnew(C2H5O2 ) * xnew(NO )  &
        + rct(22,k) * xnew(CH3COO2 ) * xnew(NO )  &
        + rct(26,k) * xnew(NO ) * xnew(SECC4H9O2 )  &
        + rct(26,k) * xnew(MEKO2 ) * xnew(NO )  &
        + rct(26,k) * xnew(ETRO2 ) * xnew(NO )  &
        + rct(26,k) * xnew(NO ) * xnew(PRRO2 )  &
        + rct(26,k) * xnew(OXYO2 ) * xnew(NO )  &
        + rct(26,k) * xnew(MALO2 ) * xnew(NO )  &
        + 0.86*rct(26,k) * xnew(ISRO2 ) * xnew(NO )  &
        + 0.95*rct(26,k) * xnew(MVKO2 ) * xnew(NO )  &
        + 0.95*rct(39,k) * xnew(MACRO2 ) * xnew(NO )  &
        + 2.5e-12 * xnew(MACRO2 ) * xnew(NO3 )  &
        + rct(41,k) * xnew(MACO3 ) * xnew(NO )  &
        + 1.9*rct(26,k) * xnew(ISNIR ) * xnew(NO )  &
        + 1.1*rct(26,k) * xnew(ISONO3 ) * xnew(NO )  &
        + 2.9e-11 * xnew(MPAN ) * xnew(OH )  &
        + rcmisc(20,k) * xnew(OP ) * xnew(NO )  &
        + rcmisc(22,k) * xnew(N2O5 )  &
        + rcmisc(25,k) * xnew(PAN )  &
        + rcmisc(29,k) * xnew(MPAN )  &
        + rcphot(IDHNO3,K) * xnew(HNO3 )  &
        + rcphot(IDNO3,K) * xnew(NO3 )  &
        + rcemis(NO2,k) 

      L =  &
          rct(2,k)* xnew(O3 )  &
        + rct(6,k)* xnew(NO3 )  &
        + rcmisc(21,k)* xnew(NO3 )  &
        + rcmisc(23,k)* xnew(OH )  &
        + rcmisc(24,k)* xnew(CH3COO2 )  &
        + rcmisc(28,k)* xnew(MACO3 )  &
        + rcphot(IDNO2,K) 

      xnew(NO2)=  max(0.0, ( xold(NO2) + dt2 * P)) /(1.0 + dt2*L )

!-> PAN 

      P =   &
          rcmisc(24,k) * xnew(CH3COO2 ) * xnew(NO2 ) 

      L =  &
          rcmisc(25,k) 

      xnew(PAN)=  max(0.0, ( xold(PAN) + dt2 * P)) /(1.0 + dt2*L )

!-> MPAN 

      P =   &
          rcmisc(28,k) * xnew(MACO3 ) * xnew(NO2 ) 

      L =  &
          2.9e-11* xnew(OH )  &
        + rcmisc(29,k) 

      xnew(MPAN)=  max(0.0, ( xold(MPAN) + dt2 * P)) /(1.0 + dt2*L )

!-> NO3 

      P =   &
          rct(2,k) * xnew(O3 ) * xnew(NO2 )  &
        + rcmisc(6,k) * xnew(OH ) * xnew(HNO3 )  &
        + rcmisc(22,k) * xnew(N2O5 ) 

      L =  &
          rct(4,k)* xnew(NO )  &
        + rct(6,k)* xnew(NO2 )  &
        + rct(16,k)* xnew(HCHO )  &
        + 2.5e-12* xnew(MACRO2 )  &
        + rct(43,k)* xnew(C5H8 )  &
        + rcmisc(21,k)* xnew(NO2 )  &
        + rcphot(IDNO3,K) 

      xnew(NO3)=  max(0.0, ( xold(NO3) + dt2 * P)) /(1.0 + dt2*L )

!-> N2O5 

      P =   &
          rcmisc(21,k) * xnew(NO2 ) * xnew(NO3 ) 

      L =  &
          rcmisc(17,k)  &
        + rcmisc(22,k) 

      xnew(N2O5)=  max(0.0, ( xold(N2O5) + dt2 * P)) /(1.0 + dt2*L )

!-> ISONO3 

      P =   &
          rct(43,k) * xnew(C5H8 ) * xnew(NO3 )  &
        + 2.0e-11 * xnew(ISONO3H ) * xnew(OH ) 

      L =  &
          rct(26,k)* xnew(NO )  &
        + rct(32,k)* xnew(HO2 ) 

      xnew(ISONO3)=  max(0.0, ( xold(ISONO3) + dt2 * P)) /(1.0 + dt2*L )

!-> HNO3 

      P =   &
          rct(16,k) * xnew(NO3 ) * xnew(HCHO )  &
        + 2.*rcmisc(17,k) * xnew(N2O5 )  &
        + rcmisc(23,k) * xnew(NO2 ) * xnew(OH ) 

      L =  &
          rcmisc(6,k)* xnew(OH )  &
        + rcmisc(19,k)  &
        + rcphot(IDHNO3,K) 

      xnew(HNO3)=  max(0.0, ( xold(HNO3) + dt2 * P)) /(1.0 + dt2*L )

!-> HONO 

      P =   &
          rcmisc(30,k) * xnew(OH ) * xnew(NO ) 

      L =  &
          rct(10,k)* xnew(OH )  &
        + 0.222*rcphot(IDNO2,K) 

      xnew(HONO)=  max(0.0, ( xold(HONO) + dt2 * P)) /(1.0 + dt2*L )

!-> CH3COO2 

      P =   &
          0.95*rct(21,k) * xnew(OH ) * xnew(CH3CHO )  &
        + rct(20,k) * xnew(CH3COO2H ) * xnew(OH )  &
        + rct(26,k) * xnew(MEKO2 ) * xnew(NO )  &
        + rct(33,k) * xnew(OH ) * xnew(MGLYOX )  &
        + 0.684*rct(26,k) * xnew(MVKO2 ) * xnew(NO )  &
        + 0.41*rct(37,k) * xnew(MACR ) * xnew(O3 )  &
        + rct(41,k) * xnew(MACO3 ) * xnew(NO )  &
        + 1.51e-11 * xnew(MACO2H ) * xnew(OH )  &
        + rcmisc(25,k) * xnew(PAN )  &
        + rcphot(IDRCOHCO,K) * xnew(MGLYOX )  &
        + rcphot(IDCH3O2H,K) * xnew(MEKO2H )  &
        + rcphot(IDCH3COX,K) * xnew(MEK ) 

      L =  &
          rct(22,k)* xnew(NO )  &
        + rct(23,k)* xnew(CH3O2 )  &
        + rct(24,k)* xnew(CH3COO2 )  &
        + rct(24,k)* xnew(CH3COO2 )  &
        + rct(25,k)* xnew(HO2 )  &
        + rcmisc(24,k)* xnew(NO2 ) 

      xnew(CH3COO2)=  max(0.0, ( xold(CH3COO2) + dt2 * P)) /(1.0 + dt2*L )

!-> MACR 

      P =   &
          0.67*rct(34,k) * xnew(C5H8 ) * xnew(O3 )  &
        + 0.32*rct(26,k) * xnew(ISRO2 ) * xnew(NO )  &
        + 0.1*rct(26,k) * xnew(ISONO3 ) * xnew(NO ) 

      L =  &
          rct(37,k)* xnew(O3 )  &
        + rct(38,k)* xnew(OH ) 

      xnew(MACR)=  max(0.0, ( xold(MACR) + dt2 * P)) /(1.0 + dt2*L )

!-> ISNI 

      P =   &
          0.14*rct(26,k) * xnew(ISRO2 ) * xnew(NO )  &
        + 0.05*rct(26,k) * xnew(MVKO2 ) * xnew(NO )  &
        + 0.05*rct(39,k) * xnew(MACRO2 ) * xnew(NO )  &
        + 0.05*rct(26,k) * xnew(ISNIR ) * xnew(NO )  &
        + 0.85*rct(26,k) * xnew(ISONO3 ) * xnew(NO ) 

      L =  &
          5.96e-11* xnew(OH ) 

      xnew(ISNI)=  max(0.0, ( xold(ISNI) + dt2 * P)) /(1.0 + dt2*L )

!-> ISNIR 

      P =   &
          5.96e-11 * xnew(ISNI ) * xnew(OH )  &
        + 3.7e-11 * xnew(ISNIRH ) * xnew(OH ) 

      L =  &
          rct(26,k)* xnew(NO )  &
        + rct(32,k)* xnew(HO2 ) 

      xnew(ISNIR)=  max(0.0, ( xold(ISNIR) + dt2 * P)) /(1.0 + dt2*L )

!-> GLYOX 

      P =   &
          rct(26,k) * xnew(MALO2 ) * xnew(NO )  &
        + rcphot(IDCH3O2H,K) * xnew(MALO2H )  &
        + rcemis(GLYOX,k) 

      L =  &
          rcmisc(16,k)* xnew(OH )  &
        + rcphot(IDHCOHCO,K) 

      xnew(GLYOX)=  max(0.0, ( xold(GLYOX) + dt2 * P)) /(1.0 + dt2*L )

!-> MGLYOX 

      P =   &
          rct(26,k) * xnew(OXYO2 ) * xnew(NO )  &
        + rct(26,k) * xnew(MALO2 ) * xnew(NO )  &
        + 0.266*rct(26,k) * xnew(MVKO2 ) * xnew(NO )  &
        + 0.59*rct(37,k) * xnew(MACR ) * xnew(O3 )  &
        + rct(40,k) * xnew(ACETOL ) * xnew(OH )  &
        + 0.82*rct(42,k) * xnew(MVK ) * xnew(O3 )  &
        + rcphot(IDCH3O2H,K) * xnew(OXYO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(MALO2H )  &
        + rcemis(MGLYOX,k) 

      L =  &
          rct(33,k)* xnew(OH )  &
        + rcphot(IDRCOHCO,K) 

      xnew(MGLYOX)=  max(0.0, ( xold(MGLYOX) + dt2 * P)) /(1.0 + dt2*L )

!-> MAL 

      P =   &
          rct(26,k) * xnew(OXYO2 ) * xnew(NO )  &
        + rcphot(IDCH3O2H,K) * xnew(OXYO2H ) 

      L =  &
          5.58e-11* xnew(OH ) 

      xnew(MAL)=  max(0.0, ( xold(MAL) + dt2 * P)) /(1.0 + dt2*L )

!-> MEK 

      P =   &
          0.65*rct(26,k) * xnew(NO ) * xnew(SECC4H9O2 )  &
        + 2.15e-11 * xnew(BURO2H ) * xnew(OH )  &
        + 0.65*rcphot(IDCH3O2H,K) * xnew(BURO2H )  &
        + rcemis(MEK,k) 

      L =  &
          rcmisc(15,k)* xnew(OH )  &
        + rcphot(IDCH3COX,K) 

      xnew(MEK)=  max(0.0, ( xold(MEK) + dt2 * P)) /(1.0 + dt2*L )

!-> MVK 

      P =   &
          0.26*rct(34,k) * xnew(C5H8 ) * xnew(O3 )  &
        + 0.42*rct(26,k) * xnew(ISRO2 ) * xnew(NO )  &
        + 0.05*rct(26,k) * xnew(ISONO3 ) * xnew(NO ) 

      L =  &
          rct(36,k)* xnew(OH )  &
        + rct(42,k)* xnew(O3 ) 

      xnew(MVK)=  max(0.0, ( xold(MVK) + dt2 * P)) /(1.0 + dt2*L )

!-> HCHO 

      P =   &
          rct(11,k) * xnew(CH3O2 ) * xnew(NO )  &
        + 2.*rct(12,k) * xnew(CH3O2 ) * xnew(CH3O2 )  &
        + rct(13,k) * xnew(CH3O2 ) * xnew(CH3O2 )  &
        + rcmisc(11,k) * xnew(OH ) * xnew(CH3OH )  &
        + 0.1*rct(14,k) * xnew(HO2 ) * xnew(CH3O2 )  &
        + 0.4*rct(15,k) * xnew(CH3O2H ) * xnew(OH )  &
        + rct(23,k) * xnew(CH3O2 ) * xnew(CH3COO2 )  &
        + 2.*rct(26,k) * xnew(ETRO2 ) * xnew(NO )  &
        + 1.14*rct(28,k) * xnew(C2H4 ) * xnew(O3 )  &
        + 0.545*rct(29,k) * xnew(O3 ) * xnew(C3H6 )  &
        + rct(26,k) * xnew(NO ) * xnew(PRRO2 )  &
        + 0.8*rct(34,k) * xnew(C5H8 ) * xnew(O3 )  &
        + 0.74*rct(26,k) * xnew(ISRO2 ) * xnew(NO )  &
        + 0.266*rct(26,k) * xnew(MVKO2 ) * xnew(NO )  &
        + 0.534*rct(37,k) * xnew(MACR ) * xnew(O3 )  &
        + 2.5e-12 * xnew(MACRO2 ) * xnew(NO3 )  &
        + rct(41,k) * xnew(MACO3 ) * xnew(NO )  &
        + 0.8*rct(42,k) * xnew(MVK ) * xnew(O3 )  &
        + 0.15*rct(26,k) * xnew(ISONO3 ) * xnew(NO )  &
        + 1.51e-11 * xnew(MACO2H ) * xnew(OH )  &
        + rcphot(IDCH3O2H,K) * xnew(CH3O2H )  &
        + 0.1*rcphot(IDHCOHCO,K) * xnew(GLYOX )  &
        + 1.56*rcphot(IDCH3O2H,K) * xnew(ETRO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(PRRO2H )  &
        + rcemis(HCHO,k) 

      L =  &
          rcmisc(12,k)* xnew(OH )  &
        + rct(16,k)* xnew(NO3 )  &
        + rcphot(IDACH2O,K)  &
        + rcphot(IDBCH2O,K) 

      xnew(HCHO)=  max(0.0, ( xold(HCHO) + dt2 * P)) /(1.0 + dt2*L )

!-> CH3CHO 

      P =   &
          rct(18,k) * xnew(C2H5O2 ) * xnew(NO )  &
        + 8.01e-12 * xnew(C2H5OOH ) * xnew(OH )  &
        + rcmisc(13,k) * xnew(OH ) * xnew(C2H5OH )  &
        + 0.35*rct(26,k) * xnew(NO ) * xnew(SECC4H9O2 )  &
        + rct(26,k) * xnew(MEKO2 ) * xnew(NO )  &
        + 1.38e-11 * xnew(ETRO2H ) * xnew(OH )  &
        + 0.545*rct(29,k) * xnew(O3 ) * xnew(C3H6 )  &
        + rct(26,k) * xnew(NO ) * xnew(PRRO2 )  &
        + 0.684*rct(26,k) * xnew(MVKO2 ) * xnew(NO )  &
        + 0.04*rct(42,k) * xnew(MVK ) * xnew(O3 )  &
        + 0.95*rct(26,k) * xnew(ISNIR ) * xnew(NO )  &
        + rcphot(IDCH3O2H,K) * xnew(C2H5OOH )  &
        + 0.22*rcphot(IDCH3O2H,K) * xnew(ETRO2H )  &
        + 0.35*rcphot(IDCH3O2H,K) * xnew(BURO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(PRRO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(MEKO2H )  &
        + rcemis(CH3CHO,k) 

      L =  &
          rct(21,k)* xnew(OH )  &
        + rcphot(IDCH3CHO,K) 

      xnew(CH3CHO)=  max(0.0, ( xold(CH3CHO) + dt2 * P)) /(1.0 + dt2*L )

!-> C2H6 

      P =   &
          rcemis(C2H6,k) 

      L =  &
          rct(17,k)* xnew(OH ) 

      xnew(C2H6)=  max(0.0, ( xold(C2H6) + dt2 * P)) /(1.0 + dt2*L )

!-> NC4H10 

      P =   &
          rcemis(NC4H10,k) 

      L =  &
          rcmisc(14,k)* xnew(OH ) 

      xnew(NC4H10)=  max(0.0, ( xold(NC4H10) + dt2 * P)) /(1.0 + dt2*L )

!-> C2H4 

      P =   &
          rcemis(C2H4,k) 

      L =  &
          rct(28,k)* xnew(O3 )  &
        + rcmisc(26,k)* xnew(OH ) 

      xnew(C2H4)=  max(0.0, ( xold(C2H4) + dt2 * P)) /(1.0 + dt2*L )

!-> C3H6 

      P =   &
          0.07*rct(34,k) * xnew(C5H8 ) * xnew(O3 )  &
        + rcemis(C3H6,k) 

      L =  &
          rct(29,k)* xnew(O3 )  &
        + rcmisc(27,k)* xnew(OH ) 

      xnew(C3H6)=  max(0.0, ( xold(C3H6) + dt2 * P)) /(1.0 + dt2*L )

!-> OXYL 

      P =   &
          rcemis(OXYL,k) 

      L =  &
          1.36e-11* xnew(OH ) 

      xnew(OXYL)=  max(0.0, ( xold(OXYL) + dt2 * P)) /(1.0 + dt2*L )

!-> C5H8 

      P =   &
          rcemis(C5H8,k)  &
        + RCBIO(BIO_ISOP,K) 

      L =  &
          rct(34,k)* xnew(O3 )  &
        + rct(35,k)* xnew(OH )  &
        + rct(43,k)* xnew(NO3 ) 

      xnew(C5H8)=  max(0.0, ( xold(C5H8) + dt2 * P)) /(1.0 + dt2*L )

!-> APINENE 

      P =   &
          RCBIO(BIO_TERP,K) 
      ! L = 0.0
 

      xnew(APINENE)=  max(0.0,  xold(APINENE) + dt2 * P)

!-> CH3O2H 

      P =   &
          0.9*rct(14,k) * xnew(HO2 ) * xnew(CH3O2 ) 

      L =  &
          rct(15,k)* xnew(OH )  &
        + 1.0e-5  &
        + rcphot(IDCH3O2H,K) 

      xnew(CH3O2H)=  max(0.0, ( xold(CH3O2H) + dt2 * P)) /(1.0 + dt2*L )

!-> C2H5OOH 

      P =   &
          rct(19,k) * xnew(C2H5O2 ) * xnew(HO2 ) 

      L =  &
          8.01e-12* xnew(OH )  &
        + rct(20,k)* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(C2H5OOH)=  max(0.0, ( xold(C2H5OOH) + dt2 * P)) /(1.0 + dt2*L )

!-> BURO2H 

      P =   &
          rct(27,k) * xnew(SECC4H9O2 ) * xnew(HO2 ) 

      L =  &
          rct(20,k)* xnew(OH )  &
        + 2.15e-11* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(BURO2H)=  max(0.0, ( xold(BURO2H) + dt2 * P)) /(1.0 + dt2*L )

!-> ETRO2H 

      P =   &
          1.2e-11 * xnew(ETRO2 ) * xnew(HO2 ) 

      L =  &
          1.38e-11* xnew(OH )  &
        + rct(20,k)* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(ETRO2H)=  max(0.0, ( xold(ETRO2H) + dt2 * P)) /(1.0 + dt2*L )

!-> PRRO2H 

      P =   &
          rct(30,k) * xnew(PRRO2 ) * xnew(HO2 ) 

      L =  &
          2.44e-11* xnew(OH )  &
        + rct(20,k)* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(PRRO2H)=  max(0.0, ( xold(PRRO2H) + dt2 * P)) /(1.0 + dt2*L )

!-> OXYO2H 

      P =   &
          rct(31,k) * xnew(OXYO2 ) * xnew(HO2 ) 

      L =  &
          4.2e-11* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(OXYO2H)=  max(0.0, ( xold(OXYO2H) + dt2 * P)) /(1.0 + dt2*L )

!-> MEKO2H 

      P =   &
          rct(27,k) * xnew(MEKO2 ) * xnew(HO2 ) 

      L =  &
          rct(20,k)* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(MEKO2H)=  max(0.0, ( xold(MEKO2H) + dt2 * P)) /(1.0 + dt2*L )

!-> MALO2H 

      P =   &
          rct(32,k) * xnew(MALO2 ) * xnew(HO2 ) 

      L =  &
          rct(20,k)* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(MALO2H)=  max(0.0, ( xold(MALO2H) + dt2 * P)) /(1.0 + dt2*L )

!-> MVKO2H 

      P =   &
          rct(27,k) * xnew(MVKO2 ) * xnew(HO2 ) 

      L =  &
          2.2e-11* xnew(OH ) 

      xnew(MVKO2H)=  max(0.0, ( xold(MVKO2H) + dt2 * P)) /(1.0 + dt2*L )

!-> MACROOH 

      P =   &
          rct(27,k) * xnew(MACRO2 ) * xnew(HO2 ) 

      L =  &
          2.82e-11* xnew(OH ) 

      xnew(MACROOH)=  max(0.0, ( xold(MACROOH) + dt2 * P)) /(1.0 + dt2*L )

!-> MACO3H 

      P =   &
          0.71*rct(44,k) * xnew(MACO3 ) * xnew(HO2 ) 

      L =  &
          1.87e-11* xnew(OH ) 

      xnew(MACO3H)=  max(0.0, ( xold(MACO3H) + dt2 * P)) /(1.0 + dt2*L )

!-> MACO2H 

      P =   &
          0.29*rct(44,k) * xnew(MACO3 ) * xnew(HO2 ) 

      L =  &
          1.51e-11* xnew(OH ) 

      xnew(MACO2H)=  max(0.0, ( xold(MACO2H) + dt2 * P)) /(1.0 + dt2*L )

!-> ISRO2H 

      P =   &
          rct(32,k) * xnew(ISRO2 ) * xnew(HO2 ) 

      L =  &
          7.5e-11* xnew(OH ) 

      xnew(ISRO2H)=  max(0.0, ( xold(ISRO2H) + dt2 * P)) /(1.0 + dt2*L )

!-> H2O2 

      P =   &
          rcmisc(7,k) * xnew(HO2 ) * xnew(HO2 )  &
        + rcmisc(8,k) * xnew(HO2 ) * xnew(HO2 )  &
        + 0.14*rct(28,k) * xnew(C2H4 ) * xnew(O3 )  &
        + 0.09*rct(29,k) * xnew(O3 ) * xnew(C3H6 )  &
        + 0.124*rct(37,k) * xnew(MACR ) * xnew(O3 ) 

      L =  &
          rct(8,k)* xnew(OH )  &
        + AQRCK(ICLRC1,K)* xnew(SO2 )  &
        + 1.0e-5  &
        + rcphot(IDH2O2,K) 

      xnew(H2O2)=  max(0.0, ( xold(H2O2) + dt2 * P)) /(1.0 + dt2*L )

!-> CH3COO2H 

      P =   &
          0.41*rct(25,k) * xnew(CH3COO2 ) * xnew(HO2 ) 

      L =  &
          rct(20,k)* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(CH3COO2H)=  max(0.0, ( xold(CH3COO2H) + dt2 * P)) /(1.0 + dt2*L )

!-> ISONO3H 

      P =   &
          rct(32,k) * xnew(ISONO3 ) * xnew(HO2 ) 

      L =  &
          2.0e-11* xnew(OH ) 

      xnew(ISONO3H)=  max(0.0, ( xold(ISONO3H) + dt2 * P)) /(1.0 + dt2*L )

!-> ISNIRH 

      P =   &
          rct(32,k) * xnew(ISNIR ) * xnew(HO2 ) 

      L =  &
          3.7e-11* xnew(OH ) 

      xnew(ISNIRH)=  max(0.0, ( xold(ISNIRH) + dt2 * P)) /(1.0 + dt2*L )

!-> CH3OH 

      P =   &
          rct(13,k) * xnew(CH3O2 ) * xnew(CH3O2 )  &
        + rcemis(CH3OH,k) 

      L =  &
          rcmisc(11,k)* xnew(OH ) 

      xnew(CH3OH)=  max(0.0, ( xold(CH3OH) + dt2 * P)) /(1.0 + dt2*L )

!-> C2H5OH 

      P =   &
          rcemis(C2H5OH,k) 

      L =  &
          rcmisc(13,k)* xnew(OH ) 

      xnew(C2H5OH)=  max(0.0, ( xold(C2H5OH) + dt2 * P)) /(1.0 + dt2*L )

!-> ACETOL 

      P =   &
          2.44e-11 * xnew(PRRO2H ) * xnew(OH )  &
        + 0.95*rct(39,k) * xnew(MACRO2 ) * xnew(NO )  &
        + 2.5e-12 * xnew(MACRO2 ) * xnew(NO3 )  &
        + 0.95*rct(26,k) * xnew(ISNIR ) * xnew(NO )  &
        + 2.9e-11 * xnew(MPAN ) * xnew(OH ) 

      L =  &
          rct(40,k)* xnew(OH ) 

      xnew(ACETOL)=  max(0.0, ( xold(ACETOL) + dt2 * P)) /(1.0 + dt2*L )

!-> H2 

      P =   &
          rcphot(IDBCH2O,K) * xnew(HCHO ) 

      L =  &
          rct(9,k)* xnew(OH ) 

      xnew(H2)=  max(0.0, ( xold(H2) + dt2 * P)) /(1.0 + dt2*L )

!-> CO 

      P =   &
          rcmisc(12,k) * xnew(OH ) * xnew(HCHO )  &
        + rct(16,k) * xnew(NO3 ) * xnew(HCHO )  &
        + 0.05*rct(21,k) * xnew(OH ) * xnew(CH3CHO )  &
        + 0.63*rct(28,k) * xnew(C2H4 ) * xnew(O3 )  &
        + 0.56*rct(29,k) * xnew(O3 ) * xnew(C3H6 )  &
        + rcmisc(16,k) * xnew(OH ) * xnew(GLYOX )  &
        + rcmisc(16,k) * xnew(OH ) * xnew(GLYOX )  &
        + rct(33,k) * xnew(OH ) * xnew(MGLYOX )  &
        + 0.05*rct(34,k) * xnew(C5H8 ) * xnew(O3 )  &
        + 0.82*rct(37,k) * xnew(MACR ) * xnew(O3 )  &
        + 0.95*rct(39,k) * xnew(MACRO2 ) * xnew(NO )  &
        + 0.05*rct(42,k) * xnew(MVK ) * xnew(O3 )  &
        + 2.9e-11 * xnew(MPAN ) * xnew(OH )  &
        + rcphot(IDACH2O,K) * xnew(HCHO )  &
        + rcphot(IDBCH2O,K) * xnew(HCHO )  &
        + rcphot(IDCH3CHO,K) * xnew(CH3CHO )  &
        + 1.9*rcphot(IDHCOHCO,K) * xnew(GLYOX )  &
        + rcphot(IDRCOHCO,K) * xnew(MGLYOX )  &
        + rcemis(CO,k) 

      L =  &
          rcmisc(10,k)* xnew(OH ) 

      xnew(CO)=  max(0.0, ( xold(CO) + dt2 * P)) /(1.0 + dt2*L )

!-> CH4 

      P =   &
          0.1*rct(29,k) * xnew(O3 ) * xnew(C3H6 ) 

      L =  &
          rcmisc(9,k)* xnew(OH ) 

      xnew(CH4)=  max(0.0, ( xold(CH4) + dt2 * P)) /(1.0 + dt2*L )

!-> SO2 

      P =   &
          rcemis(SO2,k) 

      L =  &
          2e-12*AQRCK(ICLOHSO2,K)* xnew(OH )  &
        + AQRCK(ICLRC1,K)* xnew(H2O2 )  &
        + AQRCK(ICLRC2,K)* xnew(O3 )  &
        + AQRCK(ICLRC3,K) 

      xnew(SO2)=  max(0.0, ( xold(SO2) + dt2 * P)) /(1.0 + dt2*L )
