!-> OD

  P = rcphot(IDO3_O1D,k) * xnew(O3)
  
  L = rct(2,k)  &
    + rct(3,k)  &
    + rct(4,k)
  
  xnew(OD) = (xold(OD) + dt2 * P) / (1.0 + dt2 * L)


!-> OP

  P = rct(2,k) * xnew(OD)  &
    + rct(3,k) * xnew(OD)  &
    + rcphot(IDO3_O3P,k) * xnew(O3)  &
    + rcphot(IDNO2,k) * xnew(NO2)  &
    + rcphot(IDNO3_NO2,k) * xnew(NO3)
  
  L = rct(1,k)  &
    + rct(76,k) * xnew(NO)
  
  xnew(OP) = (xold(OP) + dt2 * P) / (1.0 + dt2 * L)


!-> OH

  P = 2.* rct(4,k) * xnew(OD)  &
    + rct(8,k) * xnew(O3) * xnew(HO2)  &
    + rct(11,k) * xnew(NO) * xnew(HO2)  &
    + 4.0e-12 * xnew(HO2) * xnew(NO3)  &
    + 0.4* rct(30,k) * xnew(CH3OOH) * xnew(OH)  &
    + 8.01e-12 * xnew(C2H5OOH) * xnew(OH)  &
    + 0.44* rct(41,k) * xnew(CH3CO3) * xnew(HO2)  &
    + 2.15e-11 * xnew(BURO2H) * xnew(OH)  &
    + 1.38e-11 * xnew(ETRO2H) * xnew(OH)  &
    + 0.13* rct(48,k) * xnew(C2H4) * xnew(O3)  &
    + 0.36* rct(49,k) * xnew(O3) * xnew(C3H6)  &
    + 2.44e-11 * xnew(PRRO2H) * xnew(OH)  &
    + 0.75* rct(56,k) * xnew(O3) * xnew(C5H8)  &
    + rct(61,k) * xnew(OH) * xnew(ISRO2H)  &
    + rct(62,k) * xnew(OH) * xnew(ISRO2H)  &
    + rct(63,k) * xnew(O3) * xnew(ISON)  &
    + 0.38* rct(67,k) * xnew(O3) * xnew(MACR)  &
    + 0.38* rct(68,k) * xnew(O3) * xnew(MACR)  &
    + 3.00e-11 * xnew(OH) * xnew(MACROOH)  &
    + rcphot(IDCH3O2H,k) * xnew(ISRO2H)
  
  !-> OH cont.
    P = P +   2* rcphot(IDCH3O2H,k) *0.5 * xnew(MACROOH)  &
    + rct(75,k) * xnew(ISRO2)  &
    + 0.065* rcphot(IDNO2,k) * xnew(HPALD)  &
    + 7.61e-11 * xnew(OH) * xnew(HPALD)  &
    + 0.13* rcphot(IDNO2,k) * xnew(PACALD)  &
    + 2* rcphot(IDH2O2,k) * xnew(H2O2)  &
    + rcphot(IDNO2,k) *0.222 * xnew(HONO)  &
    + rcphot(IDHNO3,k) * xnew(HNO3)  &
    + 0.333* rcphot(IDHO2NO2,k) * xnew(HO2NO2)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(C2H5OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(ETRO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(MEKO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3CO3H)  &
    + rcphot(IDCH3O2H,k) * xnew(OXYO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(MALO2H)  &
    + 0.25* 1.0e6 * xnew(XMTO3_RO2)
  
  L = rct(7,k) * xnew(O3)  &
    + rct(13,k) * xnew(HO2)  &
    + rct(14,k) * xnew(H2O2)  &
    + rct(15,k) * xnew(H2)  &
    + rct(16,k) * xnew(HNO3)  &
    + 2.0e-11 * xnew(NO3)  &
    + rct(19,k) * xnew(HO2NO2)  &
    + rct(22,k) * xnew(HONO)  &
    + rct(23,k) * xnew(CH4)  &
    + 1.44e-13 * xnew(CO)  &
    + rct(24,k) * xnew(CO)  &
    + rct(28,k) * xnew(CH3OH)  &
    + rct(30,k) * xnew(CH3OOH)  &
    + rct(31,k) * xnew(HCHO)  &
    + rct(33,k) * xnew(C2H6)  &
    + 8.01e-12 * xnew(C2H5OOH)  &
    + rct(36,k) * xnew(C2H5OOH)  &
    + rct(37,k) * xnew(CH3CHO)  &
    + rct(42,k) * xnew(CH3CO3H)  &
    + rct(43,k) * xnew(C2H5OH)
  
  !-> OH cont.
    L = L +   rct(44,k) * xnew(NC4H10)  &
    + rct(47,k) * xnew(MEK)  &
    + rct(42,k) * xnew(MEKO2H)  &
    + rct(36,k) * xnew(BURO2H)  &
    + 2.15e-11 * xnew(BURO2H)  &
    + 1.38e-11 * xnew(ETRO2H)  &
    + rct(36,k) * xnew(ETRO2H)  &
    + 2.44e-11 * xnew(PRRO2H)  &
    + rct(36,k) * xnew(PRRO2H)  &
    + 1.36e-11 * xnew(OXYL)  &
    + 4.2e-11 * xnew(OXYO2H)  &
    + 5.58e-11 * xnew(MAL)  &
    + rct(42,k) * xnew(MALO2H)  &
    + rct(53,k) * xnew(GLYOX)  &
    + rct(54,k) * xnew(MGLYOX)  &
    + rct(55,k) * xnew(C5H8)  &
    + rct(61,k) * xnew(ISRO2H)  &
    + rct(62,k) * xnew(ISRO2H)  &
    + 9.13e-12 * xnew(IEPOX)  &
    + 3.34e-11 * xnew(ISON)
  
  !-> OH cont.
    L = L +   rct(65,k) * xnew(MACR)  &
    + rct(66,k) * xnew(MACR)  &
    + 2.90e-11 * xnew(MPAN)  &
    + 3.00e-11 * xnew(MACROOH)  &
    + rct(71,k) * xnew(ACETOL)  &
    + 1.55e-12 * xnew(NALD)  &
    + 7.61e-11 * xnew(HPALD)  &
    + rct(79,k) * xnew(NO2)  &
    + rct(81,k) * xnew(C2H4)  &
    + rct(82,k) * xnew(C3H6)  &
    + rct(83,k) * xnew(NO)  &
    + rct(79,k) * xnew(SHIPNOX)  &
    + rct(84,k) * xnew(APINENE)  &
    + rct(85,k) * xnew(BPINENE)  &
    + rct(86,k) * xnew(XTERP)  &
    + rct(88,k) * xnew(TERPOOH)  &
    + AQRCK(ICLOHSO2,K)*2.00e-12 * xnew(SO2)
  
  xnew(OH) = (xold(OH) + dt2 * P) / (1.0 + dt2 * L)


!-> HO2

  P = rct(7,k) * xnew(O3) * xnew(OH)  &
    + rct(14,k) * xnew(OH) * xnew(H2O2)  &
    + rct(15,k) * xnew(OH) * xnew(H2)  &
    + 2.0e-11 * xnew(OH) * xnew(NO3)  &
    + rct(18,k) * xnew(HO2NO2)  &
    + 1.44e-13 * xnew(OH) * xnew(CO)  &
    + rct(24,k) * xnew(OH) * xnew(CO)  &
    + rct(25,k) * xnew(CH3O2) * xnew(NO)  &
    + 2.* rct(26,k) * xnew(CH3O2) * xnew(CH3O2)  &
    + rct(28,k) * xnew(OH) * xnew(CH3OH)  &
    + rct(31,k) * xnew(OH) * xnew(HCHO)  &
    + rct(32,k) * xnew(NO3) * xnew(HCHO)  &
    + 1.2e-12 * xnew(CH3O2) * xnew(NO3)  &
    + rct(34,k) * xnew(C2H5O2) * xnew(NO)  &
    + 0.9* rct(39,k) * xnew(CH3O2) * xnew(CH3CO3)  &
    + rct(43,k) * xnew(OH) * xnew(C2H5OH)  &
    + 0.65* rct(45,k) * xnew(NO) * xnew(SECC4H9O2)  &
    + rct(45,k) * xnew(ETRO2) * xnew(NO)  &
    + 0.13* rct(48,k) * xnew(C2H4) * xnew(O3)  &
    + 0.28* rct(49,k) * xnew(O3) * xnew(C3H6)
  
  !-> HO2 cont.
    P = P +   rct(45,k) * xnew(NO) * xnew(PRRO2)  &
    + rct(45,k) * xnew(OXYO2) * xnew(NO)  &
    + rct(45,k) * xnew(MALO2) * xnew(NO)  &
    + rct(53,k) * xnew(OH) * xnew(GLYOX)  &
    + 0.75* rct(56,k) * xnew(O3) * xnew(C5H8)  &
    + rct(57,k) * xnew(NO3) * xnew(C5H8)  &
    + rct(58,k) * xnew(NO) * xnew(ISRO2)  &
    + 2* 2.00e-12 * xnew(ISRO2) * xnew(ISRO2)  &
    + 0.78* 3.34e-11 * xnew(OH) * xnew(ISON)  &
    + 0.64* rct(67,k) * xnew(O3) * xnew(MACR)  &
    + 0.64* rct(68,k) * xnew(O3) * xnew(MACR)  &
    + 1.5* rct(69,k) * xnew(NO) * xnew(MACRO2)  &
    + 1.00e-12 * xnew(MACRO2) * xnew(MACRO2)  &
    + 1.00e-12 * xnew(MACRO2) * xnew(MACRO2)  &
    + 1.17* 1.0e-12 * xnew(MACRO2) * xnew(CH3O2)  &
    + rct(71,k) * xnew(OH) * xnew(ACETOL)  &
    + rcphot(IDCH3O2H,k) * xnew(ISRO2H)  &
    + rcphot(IDCH3O2H,k) *0.64 * xnew(ISON)  &
    + rcphot(IDHCHO_H,k) *2 * xnew(MACR)  &
    + 2* rcphot(IDCH3O2H,k) *0.5 * xnew(MACROOH)
  
  !-> HO2 cont.
    P = P +   rcphot(IDMEK,k) * xnew(ACETOL)  &
    + rcphot(IDCH3CHO,k) * xnew(NALD)  &
    + rct(74,k) * xnew(ISRO2)  &
    + 0.065* rcphot(IDNO2,k) * xnew(HPALD)  &
    + 0.13* rcphot(IDNO2,k) * xnew(PACALD)  &
    + 2.0* 1.0e-12 * xnew(ISRO2) * xnew(CH3O2)  &
    + 0.667* rcphot(IDHO2NO2,k) * xnew(HO2NO2)  &
    + 2* rcphot(IDHCHO_H,k) * xnew(HCHO)  &
    + rcphot(IDCH3CHO,k) * xnew(CH3CHO)  &
    + 0.5* rcphot(IDCHOCHO,k) * xnew(GLYOX)  &
    + rcphot(IDRCOCHO,k) * xnew(MGLYOX)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(C2H5OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(ETRO2H)  &
    + 0.65* rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(PRRO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(OXYO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(MALO2H)  &
    + rct(45,k) * xnew(TERPO2) * xnew(NO)  &
    + 0.25* 1.0e6 * xnew(XMTO3_RO2)
  
  !-> HO2 cont.
    P = P +   AQRCK(ICLOHSO2,K)*2.00e-12 * xnew(OH) * xnew(SO2)
  
  L = rct(8,k) * xnew(O3)  &
    + rct(11,k) * xnew(NO)  &
    + rct(13,k) * xnew(OH)  &
    + rct(17,k) * xnew(NO2)  &
    + 4.0e-12 * xnew(NO3)  &
    + rct(20,k) * xnew(HO2)  &
    + rct(20,k) * xnew(HO2)  &
    + rct(21,k) * xnew(HO2)  &
    + rct(21,k) * xnew(HO2)  &
    + rct(29,k) * xnew(CH3O2)  &
    + rct(35,k) * xnew(C2H5O2)  &
    + rct(41,k) * xnew(CH3CO3)  &
    + rct(46,k) * xnew(SECC4H9O2)  &
    + rct(46,k) * xnew(MEKO2)  &
    + 1.2e-11 * xnew(ETRO2)  &
    + rct(50,k) * xnew(PRRO2)  &
    + rct(51,k) * xnew(OXYO2)  &
    + rct(52,k) * xnew(MALO2)  &
    + rct(60,k) * xnew(ISRO2)  &
    + rct(70,k) * xnew(MACRO2)
  
  !-> HO2 cont.
    L = L +   rct(87,k) * xnew(TERPO2)  &
    + rct(98,k)
  
  xnew(HO2) = (xold(HO2) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3O2

  P = rct(23,k) * xnew(OH) * xnew(CH4)  &
    + 0.6* rct(30,k) * xnew(CH3OOH) * xnew(OH)  &
    + 0.05* rct(37,k) * xnew(OH) * xnew(CH3CHO)  &
    + rct(38,k) * xnew(CH3CO3) * xnew(NO)  &
    + 0.9* rct(39,k) * xnew(CH3O2) * xnew(CH3CO3)  &
    + rct(40,k) * xnew(CH3CO3) * xnew(CH3CO3)  &
    + rct(40,k) * xnew(CH3CO3) * xnew(CH3CO3)  &
    + 0.44* rct(41,k) * xnew(CH3CO3) * xnew(HO2)  &
    + 0.28* rct(49,k) * xnew(O3) * xnew(C3H6)  &
    + 0.24* rct(56,k) * xnew(O3) * xnew(C5H8)  &
    + rcphot(IDCH3CHO,k) * xnew(CH3CHO)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3CO3H)  &
    + 0.08* 1.0e6 * xnew(XMTO3_RO2)
  
  L = rct(25,k) * xnew(NO)  &
    + rct(26,k) * xnew(CH3O2)  &
    + rct(26,k) * xnew(CH3O2)  &
    + rct(27,k) * xnew(CH3O2)  &
    + rct(27,k) * xnew(CH3O2)  &
    + rct(29,k) * xnew(HO2)  &
    + 1.2e-12 * xnew(NO3)  &
    + rct(39,k) * xnew(CH3CO3)  &
    + 1.0e-12 * xnew(MACRO2)  &
    + 1.0e-12 * xnew(MACRO2)  &
    + 1.0e-12 * xnew(ISRO2)  &
    + 1.0e-12 * xnew(ISRO2)
  
  xnew(CH3O2) = (xold(CH3O2) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5O2

  P = rct(33,k) * xnew(OH) * xnew(C2H6)  &
    + rct(36,k) * xnew(C2H5OOH) * xnew(OH)  &
    + 0.35* rct(45,k) * xnew(NO) * xnew(SECC4H9O2)  &
    + rcphot(IDMEK,k) * xnew(MEK)  &
    + 0.35* rcphot(IDCH3O2H,k) * xnew(BURO2H)
  
  L = rct(34,k) * xnew(NO)  &
    + rct(35,k) * xnew(HO2)
  
  xnew(C2H5O2) = (xold(C2H5O2) + dt2 * P) / (1.0 + dt2 * L)


!-> SECC4H9O2

  P = rct(44,k) * xnew(OH) * xnew(NC4H10)  &
    + rct(36,k) * xnew(BURO2H) * xnew(OH)
  
  L = rct(45,k) * xnew(NO)  &
    + rct(46,k) * xnew(HO2)
  
  xnew(SECC4H9O2) = (xold(SECC4H9O2) + dt2 * P) / (1.0 + dt2 * L)


!-> ISRO2

  P = rct(55,k) * xnew(OH) * xnew(C5H8)
  
  L = rct(58,k) * xnew(NO)  &
    + rct(59,k) * xnew(NO)  &
    + rct(60,k) * xnew(HO2)  &
    + 2.00e-12 * xnew(ISRO2)  &
    + 2.00e-12 * xnew(ISRO2)  &
    + rct(74,k)  &
    + rct(75,k)  &
    + 1.0e-12 * xnew(CH3O2)  &
    + 1.0e-12 * xnew(CH3O2)
  
  xnew(ISRO2) = (xold(ISRO2) + dt2 * P) / (1.0 + dt2 * L)


!-> ETRO2

  P = rct(36,k) * xnew(ETRO2H) * xnew(OH)  &
    + rct(81,k) * xnew(C2H4) * xnew(OH)
  
  L = rct(45,k) * xnew(NO)  &
    + 1.2e-11 * xnew(HO2)
  
  xnew(ETRO2) = (xold(ETRO2) + dt2 * P) / (1.0 + dt2 * L)


!-> PRRO2

  P = rct(36,k) * xnew(PRRO2H) * xnew(OH)  &
    + rct(82,k) * xnew(OH) * xnew(C3H6)
  
  L = rct(45,k) * xnew(NO)  &
    + rct(50,k) * xnew(HO2)
  
  xnew(PRRO2) = (xold(PRRO2) + dt2 * P) / (1.0 + dt2 * L)


!-> OXYO2

  P = 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + 4.2e-11 * xnew(OXYO2H) * xnew(OH)
  
  L = rct(45,k) * xnew(NO)  &
    + rct(51,k) * xnew(HO2)
  
  xnew(OXYO2) = (xold(OXYO2) + dt2 * P) / (1.0 + dt2 * L)


!-> MEKO2

  P = rct(47,k) * xnew(OH) * xnew(MEK)  &
    + rct(42,k) * xnew(MEKO2H) * xnew(OH)
  
  L = rct(45,k) * xnew(NO)  &
    + rct(46,k) * xnew(HO2)
  
  xnew(MEKO2) = (xold(MEKO2) + dt2 * P) / (1.0 + dt2 * L)


!-> MALO2

  P = 5.58e-11 * xnew(MAL) * xnew(OH)  &
    + rct(42,k) * xnew(MALO2H) * xnew(OH)
  
  L = rct(45,k) * xnew(NO)  &
    + rct(52,k) * xnew(HO2)
  
  xnew(MALO2) = (xold(MALO2) + dt2 * P) / (1.0 + dt2 * L)


!-> MACRO2

  P = 0.3* rct(56,k) * xnew(O3) * xnew(C5H8)  &
    + 9.13e-12 * xnew(OH) * xnew(IEPOX)  &
    + rct(65,k) * xnew(OH) * xnew(MACR)  &
    + rct(66,k) * xnew(OH) * xnew(MACR)  &
    + rct(73,k) * xnew(MPAN)  &
    + 1.0e6 * xnew(XMTO3_RO2)
  
  L = rct(69,k) * xnew(NO)  &
    + rct(69,k) * xnew(NO)  &
    + rct(70,k) * xnew(HO2)  &
    + 1.00e-12 * xnew(MACRO2)  &
    + 1.00e-12 * xnew(MACRO2)  &
    + 1.00e-12 * xnew(MACRO2)  &
    + 1.00e-12 * xnew(MACRO2)  &
    + 1.0e-12 * xnew(CH3O2)  &
    + 1.0e-12 * xnew(CH3O2)  &
    + rct(72,k) * xnew(NO2)
  
  xnew(MACRO2) = (xold(MACRO2) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3CO3

  P = 0.95* rct(37,k) * xnew(OH) * xnew(CH3CHO)  &
    + rct(42,k) * xnew(CH3CO3H) * xnew(OH)  &
    + rct(45,k) * xnew(MEKO2) * xnew(NO)  &
    + rct(54,k) * xnew(OH) * xnew(MGLYOX)  &
    + 0.3* rct(56,k) * xnew(O3) * xnew(C5H8)  &
    + 0.2* rct(67,k) * xnew(O3) * xnew(MACR)  &
    + 0.2* rct(68,k) * xnew(O3) * xnew(MACR)  &
    + 0.5* rct(69,k) * xnew(NO) * xnew(MACRO2)  &
    + 0.25* 1.0e-12 * xnew(MACRO2) * xnew(CH3O2)  &
    + rcphot(IDHCHO_H,k) *2 * xnew(MACR)  &
    + rcphot(IDMEK,k) * xnew(ACETOL)  &
    + rct(73,k) * xnew(PAN)  &
    + rcphot(IDMEK,k) * xnew(MEK)  &
    + rcphot(IDRCOCHO,k) * xnew(MGLYOX)  &
    + rcphot(IDCH3O2H,k) * xnew(MEKO2H)  &
    + 1.0e6 * xnew(XMTO3_RO2)
  
  L = rct(38,k) * xnew(NO)  &
    + rct(39,k) * xnew(CH3O2)  &
    + rct(40,k) * xnew(CH3CO3)  &
    + rct(40,k) * xnew(CH3CO3)  &
    + rct(41,k) * xnew(HO2)  &
    + rct(80,k) * xnew(NO2)
  
  xnew(CH3CO3) = (xold(CH3CO3) + dt2 * P) / (1.0 + dt2 * L)


!-> TERPO2

  P = rct(84,k) * xnew(APINENE) * xnew(OH)  &
    + rct(85,k) * xnew(BPINENE) * xnew(OH)  &
    + rct(86,k) * xnew(XTERP) * xnew(OH)  &
    + rct(88,k) * xnew(TERPOOH) * xnew(OH)  &
    + rct(92,k) * xnew(APINENE) * xnew(NO3)  &
    + 2.51e-12 * xnew(BPINENE) * xnew(NO3)  &
    + rct(93,k) * xnew(XTERP) * xnew(NO3)
  
  L = rct(45,k) * xnew(NO)  &
    + rct(87,k) * xnew(HO2)
  
  xnew(TERPO2) = (xold(TERPO2) + dt2 * P) / (1.0 + dt2 * L)


!-> XMTO3_RO2

  P = rct(89,k) * xnew(APINENE) * xnew(O3)  &
    + rct(90,k) * xnew(BPINENE) * xnew(O3)  &
    + rct(91,k) * xnew(XTERP) * xnew(O3)
  
  L = 1.0e6
  
  xnew(XMTO3_RO2) = (xold(XMTO3_RO2) + dt2 * P) / (1.0 + dt2 * L)


!-> NMVOC

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for NMVOC! xnew(NMVOC) = max(0.0, xold(NMVOC))


!-> O3

  P = rct(1,k) * xnew(OP)  &
    + 0.15* rct(41,k) * xnew(CH3CO3) * xnew(HO2)
  
  L = rct(5,k) * xnew(NO)  &
    + rct(6,k) * xnew(NO2)  &
    + rct(7,k) * xnew(OH)  &
    + rct(8,k) * xnew(HO2)  &
    + rct(48,k) * xnew(C2H4)  &
    + rct(49,k) * xnew(C3H6)  &
    + rct(56,k) * xnew(C5H8)  &
    + rct(56,k) * xnew(C5H8)  &
    + rct(56,k) * xnew(C5H8)  &
    + rct(63,k) * xnew(ISON)  &
    + rct(64,k) * xnew(ISON)  &
    + rct(67,k) * xnew(MACR)  &
    + rct(67,k) * xnew(MACR)  &
    + rct(68,k) * xnew(MACR)  &
    + rct(68,k) * xnew(MACR)  &
    + rcphot(IDO3_O1D,k)  &
    + rcphot(IDO3_O3P,k)  &
    + rct(89,k) * xnew(APINENE)  &
    + rct(90,k) * xnew(BPINENE)  &
    + rct(91,k) * xnew(XTERP)
  
  !-> O3 cont.
    L = L +   AQRCK(ICLRC2,K) * xnew(SO2)
  
  xnew(O3) = (xold(O3) + dt2 * P) / (1.0 + dt2 * L)


!-> NO

  P = rct(12,k) * xnew(NO2) * xnew(NO3)  &
    + rcphot(IDNO2,k) * xnew(NO2)  &
    + rcphot(IDNO2,k) *0.222 * xnew(HONO)  &
    + rcemis(NO,k)
  
  L = rct(5,k) * xnew(O3)  &
    + rct(9,k) * xnew(NO3)  &
    + rct(10,k) * xnew(NO)  &
    + rct(10,k) * xnew(NO)  &
    + rct(11,k) * xnew(HO2)  &
    + rct(25,k) * xnew(CH3O2)  &
    + rct(34,k) * xnew(C2H5O2)  &
    + rct(38,k) * xnew(CH3CO3)  &
    + rct(45,k) * xnew(SECC4H9O2)  &
    + rct(45,k) * xnew(MEKO2)  &
    + rct(45,k) * xnew(ETRO2)  &
    + rct(45,k) * xnew(PRRO2)  &
    + rct(45,k) * xnew(OXYO2)  &
    + rct(45,k) * xnew(MALO2)  &
    + rct(58,k) * xnew(ISRO2)  &
    + rct(59,k) * xnew(ISRO2)  &
    + rct(69,k) * xnew(MACRO2)  &
    + rct(69,k) * xnew(MACRO2)  &
    + rct(76,k) * xnew(OP)  &
    + rct(83,k) * xnew(OH)
  
  !-> NO cont.
    L = L +   rct(45,k) * xnew(TERPO2)
  
  xnew(NO) = (xold(NO) + dt2 * P) / (1.0 + dt2 * L)


!-> NO2

  P = rct(5,k) * xnew(O3) * xnew(NO)  &
    + rct(9,k) * xnew(NO) * xnew(NO3)  &
    + rct(9,k) * xnew(NO) * xnew(NO3)  &
    + rct(10,k) * xnew(NO) * xnew(NO)  &
    + rct(10,k) * xnew(NO) * xnew(NO)  &
    + rct(11,k) * xnew(NO) * xnew(HO2)  &
    + rct(12,k) * xnew(NO2) * xnew(NO3)  &
    + 2.0e-11 * xnew(OH) * xnew(NO3)  &
    + rct(18,k) * xnew(HO2NO2)  &
    + rct(19,k) * xnew(OH) * xnew(HO2NO2)  &
    + 4.0e-12 * xnew(HO2) * xnew(NO3)  &
    + rct(22,k) * xnew(OH) * xnew(HONO)  &
    + rct(25,k) * xnew(CH3O2) * xnew(NO)  &
    + 1.2e-12 * xnew(CH3O2) * xnew(NO3)  &
    + rct(34,k) * xnew(C2H5O2) * xnew(NO)  &
    + rct(38,k) * xnew(CH3CO3) * xnew(NO)  &
    + rct(45,k) * xnew(NO) * xnew(SECC4H9O2)  &
    + rct(45,k) * xnew(MEKO2) * xnew(NO)  &
    + rct(45,k) * xnew(ETRO2) * xnew(NO)  &
    + rct(45,k) * xnew(NO) * xnew(PRRO2)
  
  !-> NO2 cont.
    P = P +   rct(45,k) * xnew(OXYO2) * xnew(NO)  &
    + rct(45,k) * xnew(MALO2) * xnew(NO)  &
    + rct(58,k) * xnew(NO) * xnew(ISRO2)  &
    + 0.22* 3.34e-11 * xnew(OH) * xnew(ISON)  &
    + rct(64,k) * xnew(O3) * xnew(ISON)  &
    + 2.0* rct(69,k) * xnew(NO) * xnew(MACRO2)  &
    + 2.90e-11 * xnew(OH) * xnew(MPAN)  &
    + 1.55e-12 * xnew(OH) * xnew(NALD)  &
    + rcphot(IDCH3O2H,k) *0.64 * xnew(ISON)  &
    + rcphot(IDCH3CHO,k) * xnew(NALD)  &
    + rct(73,k) * xnew(MPAN)  &
    + rct(76,k) * xnew(OP) * xnew(NO)  &
    + rct(78,k) * xnew(N2O5)  &
    + rct(73,k) * xnew(PAN)  &
    + rcphot(IDNO3_NO2,k) * xnew(NO3)  &
    + rcphot(IDHNO3,k) * xnew(HNO3)  &
    + 0.667* rcphot(IDHO2NO2,k) * xnew(HO2NO2)  &
    + rcemis(NO2,k)  &
    + rct(45,k) * xnew(TERPO2) * xnew(NO)
  
  L = rct(6,k) * xnew(O3)  &
    + rct(12,k) * xnew(NO3)  &
    + rct(17,k) * xnew(HO2)  &
    + rct(72,k) * xnew(MACRO2)  &
    + rct(77,k) * xnew(NO3)  &
    + rct(79,k) * xnew(OH)  &
    + rct(80,k) * xnew(CH3CO3)  &
    + rcphot(IDNO2,k)
  
  xnew(NO2) = (xold(NO2) + dt2 * P) / (1.0 + dt2 * L)


!-> PAN

  P = rct(80,k) * xnew(CH3CO3) * xnew(NO2)
  
  L = rct(73,k)
  
  xnew(PAN) = (xold(PAN) + dt2 * P) / (1.0 + dt2 * L)


!-> MPAN

  P = rct(72,k) * xnew(MACRO2) * xnew(NO2)
  
  L = 2.90e-11 * xnew(OH)  &
    + rct(73,k)
  
  xnew(MPAN) = (xold(MPAN) + dt2 * P) / (1.0 + dt2 * L)


!-> NO3

  P = rct(6,k) * xnew(O3) * xnew(NO2)  &
    + rct(16,k) * xnew(OH) * xnew(HNO3)  &
    + rct(78,k) * xnew(N2O5)  &
    + 0.333* rcphot(IDHO2NO2,k) * xnew(HO2NO2)
  
  L = rct(9,k) * xnew(NO)  &
    + rct(12,k) * xnew(NO2)  &
    + 2.0e-11 * xnew(OH)  &
    + 4.0e-12 * xnew(HO2)  &
    + rct(32,k) * xnew(HCHO)  &
    + 1.2e-12 * xnew(CH3O2)  &
    + rct(57,k) * xnew(C5H8)  &
    + rct(77,k) * xnew(NO2)  &
    + rcphot(IDNO3_NO2,k)  &
    + rct(92,k) * xnew(APINENE)  &
    + 2.51e-12 * xnew(BPINENE)  &
    + rct(93,k) * xnew(XTERP)  &
    + rct(95,k)
  
  xnew(NO3) = (xold(NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> N2O5

  P = rct(77,k) * xnew(NO2) * xnew(NO3)
  
  L = rct(78,k)  &
    + rct(94,k)
  
  xnew(N2O5) = (xold(N2O5) + dt2 * P) / (1.0 + dt2 * L)


!-> HNO3

  P = rct(32,k) * xnew(NO3) * xnew(HCHO)  &
    + rct(79,k) * xnew(NO2) * xnew(OH)  &
    + rct(79,k) * xnew(SHIPNOX) * xnew(OH)  &
    + 3.2e-5 * xnew(SHIPNOX)  &
    + 2* rct(94,k) * xnew(N2O5)  &
    + rct(95,k) * xnew(NO3)
  
  L = rct(16,k) * xnew(OH)  &
    + rcphot(IDHNO3,k)  &
    + rct(96,k)  &
    + rct(97,k)
  
  xnew(HNO3) = (xold(HNO3) + dt2 * P) / (1.0 + dt2 * L)


!-> HONO

  P = rct(83,k) * xnew(OH) * xnew(NO)
  
  L = rct(22,k) * xnew(OH)  &
    + rcphot(IDNO2,k) *0.222
  
  xnew(HONO) = (xold(HONO) + dt2 * P) / (1.0 + dt2 * L)


!-> HO2NO2

  P = rct(17,k) * xnew(HO2) * xnew(NO2)
  
  L = rct(18,k)  &
    + rct(19,k) * xnew(OH)  &
    + rcphot(IDHO2NO2,k)
  
  xnew(HO2NO2) = (xold(HO2NO2) + dt2 * P) / (1.0 + dt2 * L)


!-> MACR

  P = 1.95* rct(56,k) * xnew(O3) * xnew(C5H8)  &
    + rct(58,k) * xnew(NO) * xnew(ISRO2)  &
    + 2* 2.00e-12 * xnew(ISRO2) * xnew(ISRO2)  &
    + rct(62,k) * xnew(OH) * xnew(ISRO2H)  &
    + 0.22* 3.34e-11 * xnew(OH) * xnew(ISON)  &
    + rct(64,k) * xnew(O3) * xnew(ISON)  &
    + rcphot(IDCH3O2H,k) * xnew(ISRO2H)  &
    + rcphot(IDCH3O2H,k) *0.64 * xnew(ISON)  &
    + rct(75,k) * xnew(ISRO2)  &
    + 1.0e-12 * xnew(ISRO2) * xnew(CH3O2)  &
    + rct(45,k) * xnew(TERPO2) * xnew(NO)  &
    + 0.65* 1.0e6 * xnew(XMTO3_RO2)
  
  L = rct(65,k) * xnew(OH)  &
    + rct(66,k) * xnew(OH)  &
    + rct(67,k) * xnew(O3)  &
    + rct(67,k) * xnew(O3)  &
    + rct(68,k) * xnew(O3)  &
    + rct(68,k) * xnew(O3)  &
    + rcphot(IDHCHO_H,k) *2
  
  xnew(MACR) = (xold(MACR) + dt2 * P) / (1.0 + dt2 * L)


!-> ISON

  P = rct(57,k) * xnew(NO3) * xnew(C5H8)  &
    + rct(59,k) * xnew(NO) * xnew(ISRO2)
  
  L = 3.34e-11 * xnew(OH)  &
    + rct(63,k) * xnew(O3)  &
    + rct(64,k) * xnew(O3)  &
    + rcphot(IDCH3O2H,k) *0.64
  
  xnew(ISON) = (xold(ISON) + dt2 * P) / (1.0 + dt2 * L)


!-> GLYOX

  P = rct(45,k) * xnew(MALO2) * xnew(NO)  &
    + rcphot(IDCH3O2H,k) * xnew(MALO2H)  &
    + rcemis(GLYOX,k)
  
  L = rct(53,k) * xnew(OH)  &
    + rcphot(IDCHOCHO,k)
  
  xnew(GLYOX) = (xold(GLYOX) + dt2 * P) / (1.0 + dt2 * L)


!-> MGLYOX

  P = rct(45,k) * xnew(OXYO2) * xnew(NO)  &
    + rct(45,k) * xnew(MALO2) * xnew(NO)  &
    + 1.8* rct(67,k) * xnew(O3) * xnew(MACR)  &
    + 1.8* rct(68,k) * xnew(O3) * xnew(MACR)  &
    + rct(69,k) * xnew(NO) * xnew(MACRO2)  &
    + 2.0* 1.00e-12 * xnew(MACRO2) * xnew(MACRO2)  &
    + 1.0e-12 * xnew(MACRO2) * xnew(CH3O2)  &
    + rct(71,k) * xnew(OH) * xnew(ACETOL)  &
    + rcphot(IDCH3O2H,k) *0.5 * xnew(MACROOH)  &
    + 7.61e-11 * xnew(OH) * xnew(HPALD)  &
    + 0.13* rcphot(IDNO2,k) * xnew(PACALD)  &
    + 0.5* 1.0e-12 * xnew(ISRO2) * xnew(CH3O2)  &
    + rcphot(IDCH3O2H,k) * xnew(OXYO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(MALO2H)  &
    + rcemis(MGLYOX,k)
  
  L = rct(54,k) * xnew(OH)  &
    + rcphot(IDRCOCHO,k)
  
  xnew(MGLYOX) = (xold(MGLYOX) + dt2 * P) / (1.0 + dt2 * L)


!-> MAL

  P = rct(45,k) * xnew(OXYO2) * xnew(NO)  &
    + rcphot(IDCH3O2H,k) * xnew(OXYO2H)
  
  L = 5.58e-11 * xnew(OH)
  
  xnew(MAL) = (xold(MAL) + dt2 * P) / (1.0 + dt2 * L)


!-> MEK

  P = 0.65* rct(45,k) * xnew(NO) * xnew(SECC4H9O2)  &
    + 2.15e-11 * xnew(BURO2H) * xnew(OH)  &
    + 0.65* rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + rcemis(MEK,k)
  
  L = rct(47,k) * xnew(OH)  &
    + rcphot(IDMEK,k)
  
  xnew(MEK) = (xold(MEK) + dt2 * P) / (1.0 + dt2 * L)


!-> HCHO

  P = rct(25,k) * xnew(CH3O2) * xnew(NO)  &
    + 2.* rct(26,k) * xnew(CH3O2) * xnew(CH3O2)  &
    + rct(27,k) * xnew(CH3O2) * xnew(CH3O2)  &
    + rct(28,k) * xnew(OH) * xnew(CH3OH)  &
    + 0.1* rct(29,k) * xnew(HO2) * xnew(CH3O2)  &
    + 0.4* rct(30,k) * xnew(CH3OOH) * xnew(OH)  &
    + 1.2e-12 * xnew(CH3O2) * xnew(NO3)  &
    + rct(39,k) * xnew(CH3O2) * xnew(CH3CO3)  &
    + 2.* rct(45,k) * xnew(ETRO2) * xnew(NO)  &
    + 1.14* rct(48,k) * xnew(C2H4) * xnew(O3)  &
    + 0.545* rct(49,k) * xnew(O3) * xnew(C3H6)  &
    + rct(45,k) * xnew(NO) * xnew(PRRO2)  &
    + 1.740* rct(56,k) * xnew(O3) * xnew(C5H8)  &
    + rct(58,k) * xnew(NO) * xnew(ISRO2)  &
    + 2* 2.00e-12 * xnew(ISRO2) * xnew(ISRO2)  &
    + 0.22* 3.34e-11 * xnew(OH) * xnew(ISON)  &
    + rct(64,k) * xnew(O3) * xnew(ISON)  &
    + 1.5* rct(69,k) * xnew(NO) * xnew(MACRO2)  &
    + 1.00e-12 * xnew(MACRO2) * xnew(MACRO2)  &
    + 2.75* 1.0e-12 * xnew(MACRO2) * xnew(CH3O2)
  
  !-> HCHO cont.
    P = P +   1.55e-12 * xnew(OH) * xnew(NALD)  &
    + rcphot(IDCH3O2H,k) * xnew(ISRO2H)  &
    + rcphot(IDCH3O2H,k) *0.64 * xnew(ISON)  &
    + rcphot(IDHCHO_H,k) *2 * xnew(MACR)  &
    + rcphot(IDCH3O2H,k) *0.5 * xnew(MACROOH)  &
    + rcphot(IDMEK,k) * xnew(ACETOL)  &
    + rcphot(IDCH3CHO,k) * xnew(NALD)  &
    + rct(75,k) * xnew(ISRO2)  &
    + 7.61e-11 * xnew(OH) * xnew(HPALD)  &
    + 2.5* 1.0e-12 * xnew(ISRO2) * xnew(CH3O2)  &
    + 0.1* rcphot(IDCHOCHO,k) * xnew(GLYOX)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3OOH)  &
    + 1.56* rcphot(IDCH3O2H,k) * xnew(ETRO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(PRRO2H)  &
    + rcemis(HCHO,k)  &
    + 0.58* 1.0e6 * xnew(XMTO3_RO2)
  
  L = rct(31,k) * xnew(OH)  &
    + rct(32,k) * xnew(NO3)  &
    + rcphot(IDHCHO_H,k)  &
    + rcphot(IDHCHO_H2,k)
  
  xnew(HCHO) = (xold(HCHO) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3CHO

  P = rct(34,k) * xnew(C2H5O2) * xnew(NO)  &
    + 8.01e-12 * xnew(C2H5OOH) * xnew(OH)  &
    + rct(43,k) * xnew(OH) * xnew(C2H5OH)  &
    + 0.35* rct(45,k) * xnew(NO) * xnew(SECC4H9O2)  &
    + rct(45,k) * xnew(MEKO2) * xnew(NO)  &
    + 1.38e-11 * xnew(ETRO2H) * xnew(OH)  &
    + 0.545* rct(49,k) * xnew(O3) * xnew(C3H6)  &
    + rct(45,k) * xnew(NO) * xnew(PRRO2)  &
    + rcphot(IDCH3O2H,k) * xnew(C2H5OOH)  &
    + 0.22* rcphot(IDCH3O2H,k) * xnew(ETRO2H)  &
    + 0.35* rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(PRRO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(MEKO2H)  &
    + rcemis(CH3CHO,k)
  
  L = rct(37,k) * xnew(OH)  &
    + rcphot(IDCH3CHO,k)
  
  xnew(CH3CHO) = (xold(CH3CHO) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H6

  P = rcemis(C2H6,k)
  
  L = rct(33,k) * xnew(OH)
  
  xnew(C2H6) = (xold(C2H6) + dt2 * P) / (1.0 + dt2 * L)


!-> NC4H10

  P = rcemis(NC4H10,k)
  
  L = rct(44,k) * xnew(OH)
  
  xnew(NC4H10) = (xold(NC4H10) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H4

  P = rcemis(C2H4,k)
  
  L = rct(48,k) * xnew(O3)  &
    + rct(81,k) * xnew(OH)
  
  xnew(C2H4) = (xold(C2H4) + dt2 * P) / (1.0 + dt2 * L)


!-> C3H6

  P = rcemis(C3H6,k)
  
  L = rct(49,k) * xnew(O3)  &
    + rct(82,k) * xnew(OH)
  
  xnew(C3H6) = (xold(C3H6) + dt2 * P) / (1.0 + dt2 * L)


!-> OXYL

  P = rcemis(OXYL,k)
  
  L = 1.36e-11 * xnew(OH)
  
  xnew(OXYL) = (xold(OXYL) + dt2 * P) / (1.0 + dt2 * L)


!-> C5H8

  P = rcemis(C5H8,k)  &
    + RCBIO(1,K)
  
  L = rct(55,k) * xnew(OH)  &
    + rct(56,k) * xnew(O3)  &
    + rct(56,k) * xnew(O3)  &
    + rct(56,k) * xnew(O3)  &
    + rct(57,k) * xnew(NO3)
  
  xnew(C5H8) = (xold(C5H8) + dt2 * P) / (1.0 + dt2 * L)


!-> APINENE

  P = 0.45* RCBIO(2,K)
  
  L = rct(84,k) * xnew(OH)  &
    + rct(89,k) * xnew(O3)  &
    + rct(92,k) * xnew(NO3)
  
  xnew(APINENE) = (xold(APINENE) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3OOH

  P = 0.9* rct(29,k) * xnew(HO2) * xnew(CH3O2)
  
  L = rct(30,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)  &
    + 1.0e-5
  
  xnew(CH3OOH) = (xold(CH3OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5OOH

  P = rct(35,k) * xnew(C2H5O2) * xnew(HO2)
  
  L = 8.01e-12 * xnew(OH)  &
    + rct(36,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(C2H5OOH) = (xold(C2H5OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> BURO2H

  P = rct(46,k) * xnew(SECC4H9O2) * xnew(HO2)
  
  L = rct(36,k) * xnew(OH)  &
    + 2.15e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(BURO2H) = (xold(BURO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> ETRO2H

  P = 1.2e-11 * xnew(ETRO2) * xnew(HO2)
  
  L = 1.38e-11 * xnew(OH)  &
    + rct(36,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(ETRO2H) = (xold(ETRO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> PRRO2H

  P = rct(50,k) * xnew(PRRO2) * xnew(HO2)
  
  L = 2.44e-11 * xnew(OH)  &
    + rct(36,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(PRRO2H) = (xold(PRRO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> OXYO2H

  P = 0.96* rct(51,k) * xnew(OXYO2) * xnew(HO2)
  
  L = 4.2e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(OXYO2H) = (xold(OXYO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> MEKO2H

  P = rct(46,k) * xnew(MEKO2) * xnew(HO2)
  
  L = rct(42,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(MEKO2H) = (xold(MEKO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> MALO2H

  P = rct(52,k) * xnew(MALO2) * xnew(HO2)
  
  L = rct(42,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(MALO2H) = (xold(MALO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> MACROOH

  P = rct(70,k) * xnew(HO2) * xnew(MACRO2)
  
  L = 3.00e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k) *0.5  &
    + rcphot(IDCH3O2H,k) *0.5
  
  xnew(MACROOH) = (xold(MACROOH) + dt2 * P) / (1.0 + dt2 * L)


!-> ISRO2H

  P = rct(60,k) * xnew(HO2) * xnew(ISRO2)
  
  L = rct(61,k) * xnew(OH)  &
    + rct(62,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(ISRO2H) = (xold(ISRO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> H2O2

  P = rct(20,k) * xnew(HO2) * xnew(HO2)  &
    + rct(21,k) * xnew(HO2) * xnew(HO2)  &
    + 0.14* rct(48,k) * xnew(C2H4) * xnew(O3)  &
    + 0.09* rct(49,k) * xnew(O3) * xnew(C3H6)  &
    + 0.27* rct(56,k) * xnew(O3) * xnew(C5H8)  &
    + 0.09* 1.0e6 * xnew(XMTO3_RO2)  &
    + 0.5* rct(98,k) * xnew(HO2)
  
  L = rct(14,k) * xnew(OH)  &
    + rcphot(IDH2O2,k)  &
    + AQRCK(ICLRC1,K) * xnew(SO2)  &
    + 1.0e-5
  
  xnew(H2O2) = (xold(H2O2) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3CO3H

  P = 0.41* rct(41,k) * xnew(CH3CO3) * xnew(HO2)
  
  L = rct(42,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(CH3CO3H) = (xold(CH3CO3H) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3OH

  P = rct(27,k) * xnew(CH3O2) * xnew(CH3O2)  &
    + 0.5* 1.0e-12 * xnew(ISRO2) * xnew(CH3O2)  &
    + rcemis(CH3OH,k)
  
  L = rct(28,k) * xnew(OH)
  
  xnew(CH3OH) = (xold(CH3OH) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5OH

  P = rcemis(C2H5OH,k)
  
  L = rct(43,k) * xnew(OH)
  
  xnew(C2H5OH) = (xold(C2H5OH) + dt2 * P) / (1.0 + dt2 * L)


!-> ACETOL

  P = 2.44e-11 * xnew(PRRO2H) * xnew(OH)  &
    + 0.78* 3.34e-11 * xnew(OH) * xnew(ISON)  &
    + 0.5* rct(69,k) * xnew(NO) * xnew(MACRO2)  &
    + 2.0* 1.00e-12 * xnew(MACRO2) * xnew(MACRO2)  &
    + 0.75* 1.0e-12 * xnew(MACRO2) * xnew(CH3O2)  &
    + 2.90e-11 * xnew(OH) * xnew(MPAN)  &
    + 3.00e-11 * xnew(OH) * xnew(MACROOH)  &
    + rcphot(IDCH3O2H,k) *0.5 * xnew(MACROOH)  &
    + 0.5* 1.0e-12 * xnew(ISRO2) * xnew(CH3O2)
  
  L = rct(71,k) * xnew(OH)  &
    + rcphot(IDMEK,k)
  
  xnew(ACETOL) = (xold(ACETOL) + dt2 * P) / (1.0 + dt2 * L)


!-> NALD

  P = 0.78* 3.34e-11 * xnew(OH) * xnew(ISON)  &
    + rct(63,k) * xnew(O3) * xnew(ISON)
  
  L = 1.55e-12 * xnew(OH)  &
    + rcphot(IDCH3CHO,k)
  
  xnew(NALD) = (xold(NALD) + dt2 * P) / (1.0 + dt2 * L)


!-> HPALD

  P = rct(74,k) * xnew(ISRO2)
  
  L = 0.065* rcphot(IDNO2,k)  &
    + 7.61e-11 * xnew(OH)
  
  xnew(HPALD) = (xold(HPALD) + dt2 * P) / (1.0 + dt2 * L)


!-> PACALD

  P = 0.065* rcphot(IDNO2,k) * xnew(HPALD)
  
  L = 0.13* rcphot(IDNO2,k)
  
  xnew(PACALD) = (xold(PACALD) + dt2 * P) / (1.0 + dt2 * L)


!-> IEPOX

  P = rct(61,k) * xnew(OH) * xnew(ISRO2H)
  
  L = 9.13e-12 * xnew(OH)
  
  xnew(IEPOX) = (xold(IEPOX) + dt2 * P) / (1.0 + dt2 * L)


!-> H2

  P = rcphot(IDHCHO_H2,k) * xnew(HCHO)
  
  L = rct(15,k) * xnew(OH)
  
  xnew(H2) = (xold(H2) + dt2 * P) / (1.0 + dt2 * L)


!-> CO

  P = rct(31,k) * xnew(OH) * xnew(HCHO)  &
    + rct(32,k) * xnew(NO3) * xnew(HCHO)  &
    + 0.05* rct(37,k) * xnew(OH) * xnew(CH3CHO)  &
    + 0.63* rct(48,k) * xnew(C2H4) * xnew(O3)  &
    + 0.56* rct(49,k) * xnew(O3) * xnew(C3H6)  &
    + rct(53,k) * xnew(OH) * xnew(GLYOX)  &
    + rct(53,k) * xnew(OH) * xnew(GLYOX)  &
    + rct(54,k) * xnew(OH) * xnew(MGLYOX)  &
    + 0.42* rct(56,k) * xnew(O3) * xnew(C5H8)  &
    + 0.44* rct(67,k) * xnew(O3) * xnew(MACR)  &
    + 0.44* rct(68,k) * xnew(O3) * xnew(MACR)  &
    + 0.5* rct(69,k) * xnew(NO) * xnew(MACRO2)  &
    + 1.00e-12 * xnew(MACRO2) * xnew(MACRO2)  &
    + 0.25* 1.0e-12 * xnew(MACRO2) * xnew(CH3O2)  &
    + 3.00e-11 * xnew(OH) * xnew(MACROOH)  &
    + 1.55e-12 * xnew(OH) * xnew(NALD)  &
    + rcphot(IDHCHO_H,k) *2 * xnew(MACR)  &
    + rcphot(IDCH3O2H,k) *0.5 * xnew(MACROOH)  &
    + rcphot(IDCH3CHO,k) * xnew(NALD)  &
    + 7.61e-11 * xnew(OH) * xnew(HPALD)
  
  !-> CO cont.
    P = P +   0.13* rcphot(IDNO2,k) * xnew(PACALD)  &
    + rcphot(IDHCHO_H,k) * xnew(HCHO)  &
    + rcphot(IDHCHO_H2,k) * xnew(HCHO)  &
    + rcphot(IDCH3CHO,k) * xnew(CH3CHO)  &
    + 1.9* rcphot(IDCHOCHO,k) * xnew(GLYOX)  &
    + rcphot(IDRCOCHO,k) * xnew(MGLYOX)  &
    + rcemis(CO,k)  &
    + 0.14* 1.0e6 * xnew(XMTO3_RO2)
  
  L = 1.44e-13 * xnew(OH)  &
    + rct(24,k) * xnew(OH)
  
  xnew(CO) = (xold(CO) + dt2 * P) / (1.0 + dt2 * L)


!-> CH4

  P = 0.1* rct(49,k) * xnew(O3) * xnew(C3H6)
  
  L = rct(23,k) * xnew(OH)
  
  xnew(CH4) = (xold(CH4) + dt2 * P) / (1.0 + dt2 * L)


!-> SO2

  P = rcemis(SO2,k)
  
  L = AQRCK(ICLOHSO2,K)*2.00e-12 * xnew(OH)  &
    + AQRCK(ICLRC1,K) * xnew(H2O2)  &
    + AQRCK(ICLRC2,K) * xnew(O3)  &
    + AQRCK(ICLRC3,K)
  
  xnew(SO2) = (xold(SO2) + dt2 * P) / (1.0 + dt2 * L)


!-> SHIPNOX

  P = rcemis(SHIPNOX,k)
  
  L = rct(79,k) * xnew(OH)  &
    + 3.2e-5
  
  xnew(SHIPNOX) = (xold(SHIPNOX) + dt2 * P) / (1.0 + dt2 * L)


!-> MVK

  P = rct(45,k) * xnew(TERPO2) * xnew(NO)
  
  ! L = 0.0
  
  xnew(MVK) = xold(MVK) + dt2 * P


!-> BPINENE

  P = 0.2* RCBIO(2,K)
  
  L = rct(85,k) * xnew(OH)  &
    + rct(90,k) * xnew(O3)  &
    + 2.51e-12 * xnew(NO3)
  
  xnew(BPINENE) = (xold(BPINENE) + dt2 * P) / (1.0 + dt2 * L)


!-> XTERP

  P = 0.35* RCBIO(2,K)
  
  L = rct(86,k) * xnew(OH)  &
    + rct(91,k) * xnew(O3)  &
    + rct(93,k) * xnew(NO3)
  
  xnew(XTERP) = (xold(XTERP) + dt2 * P) / (1.0 + dt2 * L)


!-> SQT_SOA_NV

  P = 0.01534* RCBIO(2,K)
  
  ! L = 0.0
  
  xnew(SQT_SOA_NV) = xold(SQT_SOA_NV) + dt2 * P


!-> TERPOOH

  P = rct(87,k) * xnew(TERPO2) * xnew(HO2)
  
  L = rct(88,k) * xnew(OH)
  
  xnew(TERPOOH) = (xold(TERPOOH) + dt2 * P) / (1.0 + dt2 * L)


!-> TERPPeroxy

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for TERPPeroxy! xnew(TERPPeroxy) = max(0.0, xold(TERPPeroxy))


!-> VBS_TEST

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for VBS_TEST! xnew(VBS_TEST) = max(0.0, xold(VBS_TEST))


!-> Dust_ROAD_f

  P = rcemis(DUST_ROAD_F,K)
  
  ! L = 0.0
  
  xnew(Dust_ROAD_f) = xold(Dust_ROAD_f) + dt2 * P


!-> Dust_ROAD_c

  P = rcemis(DUST_ROAD_C,K)
  
  ! L = 0.0
  
  xnew(Dust_ROAD_c) = xold(Dust_ROAD_c) + dt2 * P


!-> Dust_WB_f

  P = rcemis(DUST_WB_F,K)
  
  ! L = 0.0
  
  xnew(Dust_WB_f) = xold(Dust_WB_f) + dt2 * P


!-> Dust_WB_c

  P = rcemis(DUST_WB_C,K)
  
  ! L = 0.0
  
  xnew(Dust_WB_c) = xold(Dust_WB_c) + dt2 * P


!-> Dust_SAH_f

  P = rcemis(DUST_SAH_F,K)
  
  ! L = 0.0
  
  xnew(Dust_SAH_f) = xold(Dust_SAH_f) + dt2 * P


!-> Dust_SAH_c

  P = rcemis(DUST_SAH_C,K)
  
  ! L = 0.0
  
  xnew(Dust_SAH_c) = xold(Dust_SAH_c) + dt2 * P


