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
    + 0.873* rcphot(IDNO3,k) * xnew(NO3)
  
  L = rct(1,k)  &
    + rct(67,k) * xnew(NO)
  
  xnew(OP) = (xold(OP) + dt2 * P) / (1.0 + dt2 * L)


!-> OH

  P = 2* rct(4,k) * xnew(OD)  &
    + rct(8,k) * xnew(O3) * xnew(HO2)  &
    + rct(10,k) * xnew(NO) * xnew(HO2)  &
    + 4.00e-12 * xnew(HO2) * xnew(NO3)  &
    + 0.4* rct(29,k) * xnew(CH3OOH) * xnew(OH)  &
    + 0.69* 1.17e-11 * xnew(C2H5OOH) * xnew(OH)  &
    + 0.454* rct(38,k) * xnew(CH3CO3) * xnew(HO2)  &
    + 0.847* 2.4e-11 * xnew(BURO2H) * xnew(OH)  &
    + 0.826* 2.1e-11 * xnew(MEKO2H) * xnew(OH)  &
    + 0.79* 1.7475e-11 * xnew(ETRO2H) * xnew(OH)  &
    + 0.17* rct(43,k) * xnew(C2H4) * xnew(O3)  &
    + 0.3646* rct(44,k) * xnew(O3) * xnew(C3H6)  &
    + 0.87* 2.81e-11 * xnew(PRRO2H) * xnew(OH)  &
    + 0.97* 1.25e-10 * xnew(C5DICAROOH) * xnew(OH)  &
    + 0.584* rct(53,k) * xnew(OH) * xnew(GLYOX)  &
    + 0.25* rct(56,k) * xnew(O3) * xnew(C5H8)  &
    + 8.94e-11 * xnew(OH) * xnew(ISRO2H)  &
    + 0.67* 9.20e-17 * xnew(O3) * xnew(ISON)  &
    + 0.19* rct(60,k) * xnew(O3) * xnew(MACR)  &
    + 0.19* rct(61,k) * xnew(O3) * xnew(MACR)
  
  !-> OH cont.
    P = P +   3.00e-11 * xnew(OH) * xnew(MACROOH)  &
    + rcphot(IDCH3O2H,k) * xnew(ISRO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(MACROOH)  &
    + rct(66,k) * xnew(ISRO2)  &
    + rcphot(IDHCHO_H2,k) *12. * xnew(HPALD)  &
    + 7.61e-11 * xnew(OH) * xnew(HPALD)  &
    + rcphot(IDHCHO_H2,k) *24. * xnew(PACALD)  &
    + 2* rcphot(IDH2O2,k) * xnew(H2O2)  &
    + rcphot(IDHONO,k) * xnew(HONO)  &
    + rcphot(IDHNO3,k) * xnew(HNO3)  &
    + 0.333* rcphot(IDHO2NO2,k) * xnew(HO2NO2)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(C2H5OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(ETRO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(MEKO2H)  &
    + rcphot(IDMEK,k) * xnew(MEKO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3CO3H)  &
    + rcphot(IDCH3O2H,k) * xnew(OXYO2H)  &
    + rcphot(IDHCHO_H,k) *0.9 * xnew(C5DICAROOH)
  
  !-> OH cont.
    P = P +   0.36* rct(99,k) * xnew(MVK) * xnew(O3)  &
    + 0.7* rct(76,k) * xnew(APINENE) * xnew(O3)
  
  L = rct(7,k) * xnew(O3)  &
    + rct(12,k) * xnew(HO2)  &
    + rct(13,k) * xnew(H2O2)  &
    + rct(14,k) * xnew(H2)  &
    + rct(15,k) * xnew(HNO3)  &
    + rct(18,k) * xnew(HONO)  &
    + 2.00e-11 * xnew(NO3)  &
    + rct(21,k) * xnew(HO2NO2)  &
    + rct(22,k) * xnew(CH4)  &
    + rct(23,k) * xnew(CO)  &
    + rct(28,k) * xnew(CH3OH)  &
    + rct(29,k) * xnew(CH3OOH)  &
    + rct(30,k) * xnew(HCHO)  &
    + rct(32,k) * xnew(C2H6)  &
    + 1.17e-11 * xnew(C2H5OOH)  &
    + rct(35,k) * xnew(CH3CHO)  &
    + 1.1e-11 * xnew(CH3CO3H)  &
    + 3.22e-12 * xnew(C2H5OH)  &
    + rct(39,k) * xnew(NC4H10)  &
    + 2.4e-11 * xnew(BURO2H)
  
  !-> OH cont.
    L = L +   8.6e-13 * xnew(SC4H9NO3)  &
    + 1.1e-12 * xnew(MEK)  &
    + 2.1e-11 * xnew(MEKO2H)  &
    + 2.3e-13 * xnew(BIACET)  &
    + 1.7475e-11 * xnew(ETRO2H)  &
    + 2.81e-11 * xnew(PRRO2H)  &
    + rct(47,k) * xnew(BENZENE)  &
    + rct(48,k) * xnew(TOLUENE)  &
    + 1.36e-11 * xnew(OXYL)  &
    + 7.62e-11 * xnew(OXYO2H)  &
    + 6.2e-11 * xnew(C5DICARB)  &
    + 1.25e-10 * xnew(C5DICAROOH)  &
    + 7.48e-11 * xnew(C5134CO2OH)  &
    + 1.72e-11 * xnew(C54CO)  &
    + rct(53,k) * xnew(GLYOX)  &
    + rct(54,k) * xnew(MGLYOX)  &
    + rct(55,k) * xnew(C5H8)  &
    + 8.94e-11 * xnew(ISRO2H)  &
    + 9.13e-12 * xnew(IEPOX)  &
    + 3.34e-11 * xnew(ISON)
  
  !-> OH cont.
    L = L +   rct(58,k) * xnew(MACR)  &
    + rct(59,k) * xnew(MACR)  &
    + 2.90e-11 * xnew(MPAN)  &
    + 3.00e-11 * xnew(MACROOH)  &
    + rct(62,k) * xnew(ACETOL)  &
    + 1.55e-12 * xnew(NALD)  &
    + 7.61e-11 * xnew(HPALD)  &
    + rct(70,k) * xnew(NO2)  &
    + rct(72,k) * xnew(C2H4)  &
    + rct(73,k) * xnew(C3H6)  &
    + rct(74,k) * xnew(NO)  &
    + rct(95,k) * xnew(APINENE)  &
    + rct(97,k) * xnew(TERPOOH)  &
    + rct(98,k) * xnew(MVK)  &
    + AQRCK(ICLOHSO2,K)*2.00e-12 * xnew(SO2)  &
    + rct(70,k) * xnew(shipNOx)
  
  xnew(OH) = (xold(OH) + dt2 * P) / (1.0 + dt2 * L)


!-> HO2

  P = rct(7,k) * xnew(O3) * xnew(OH)  &
    + rct(13,k) * xnew(OH) * xnew(H2O2)  &
    + rct(14,k) * xnew(OH) * xnew(H2)  &
    + 2.00e-11 * xnew(OH) * xnew(NO3)  &
    + rct(20,k) * xnew(HO2NO2)  &
    + rct(23,k) * xnew(OH) * xnew(CO)  &
    + 1.2e-12 * xnew(CH3O2) * xnew(NO3)  &
    + rct(24,k) * xnew(CH3O2) * xnew(NO)  &
    + rct(25,k) * xnew(RO2POOL) * xnew(CH3O2)  &
    + rct(28,k) * xnew(CH3OH) * xnew(OH)  &
    + rct(30,k) * xnew(OH) * xnew(HCHO)  &
    + rct(31,k) * xnew(NO3) * xnew(HCHO)  &
    + rct(33,k) * xnew(C2H5O2) * xnew(NO)  &
    + 2.3e-12 * xnew(C2H5O2) * xnew(NO3)  &
    + 0.95* 3.22e-12 * xnew(OH) * xnew(C2H5OH)  &
    + 0.654* rct(41,k) * xnew(NO) * xnew(C4H9O2)  &
    + 0.716* 2.3e-12 * xnew(NO3) * xnew(C4H9O2)  &
    + 0.652* rcphot(IDMEK,k) *0.78 * xnew(SC4H9NO3)  &
    + rct(41,k) * xnew(ETRO2) * xnew(NO)  &
    + 0.27* rct(43,k) * xnew(C2H4) * xnew(O3)
  
  !-> HO2 cont.
    P = P +   0.038* rct(44,k) * xnew(O3) * xnew(C3H6)  &
    + rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 2.3e-12 * xnew(NO3) * xnew(PRRO2)  &
    + rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + 2.3e-12 * xnew(OXYO2) * xnew(NO3)  &
    + rcphot(IDNO2,k) *0.2 * xnew(C5DICARB)  &
    + rct(41,k) * xnew(C5DICARBO2) * xnew(NO)  &
    + rct(50,k) * xnew(C5DICARBO2) * xnew(NO3)  &
    + 7.48e-11 * xnew(C5134CO2OH) * xnew(OH)  &
    + rcphot(IDRCOCHO,k) * xnew(C54CO)  &
    + rcphot(IDCH3COY,k) *2 * xnew(C54CO)  &
    + 0.416* rct(53,k) * xnew(OH) * xnew(GLYOX)  &
    + 0.25* rct(56,k) * xnew(O3) * xnew(C5H8)  &
    + rct(57,k) * xnew(NO3) * xnew(C5H8)  &
    + 0.9* rct(41,k) * xnew(NO) * xnew(ISRO2)  &
    + 4.00e-12 * xnew(RO2POOL) * xnew(ISRO2)  &
    + 0.78* 3.34e-11 * xnew(OH) * xnew(ISON)  &
    + 0.32* rct(60,k) * xnew(O3) * xnew(MACR)  &
    + 0.32* rct(61,k) * xnew(O3) * xnew(MACR)  &
    + 0.75* rct(41,k) * xnew(NO) * xnew(MACRO2)
  
  !-> HO2 cont.
    P = P +   rct(62,k) * xnew(OH) * xnew(ACETOL)  &
    + rcphot(IDCH3O2H,k) * xnew(ISRO2H)  &
    + rcphot(IDMEK,k) *0.78 * xnew(ISON)  &
    + rcphot(IDHCHO_H2,k) *0.0474 * xnew(MACR)  &
    + rcphot(IDCH3O2H,k) * xnew(MACROOH)  &
    + rcphot(IDMEK,k) * xnew(ACETOL)  &
    + rcphot(IDCH3CHO,k) * xnew(NALD)  &
    + rct(65,k) * xnew(ISRO2)  &
    + rcphot(IDHCHO_H2,k) *12. * xnew(HPALD)  &
    + rcphot(IDHCHO_H2,k) *24. * xnew(PACALD)  &
    + 0.667* rcphot(IDHO2NO2,k) * xnew(HO2NO2)  &
    + 2* rcphot(IDHCHO_H,k) * xnew(HCHO)  &
    + rcphot(IDCH3CHO,k) * xnew(CH3CHO)  &
    + 0.5* rcphot(IDCHOCHO,k) * xnew(GLYOX)  &
    + rcphot(IDRCOCHO,k) * xnew(MGLYOX)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(C2H5OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(ETRO2H)  &
    + 0.716* rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(PRRO2H)
  
  !-> HO2 cont.
    P = P +   rcphot(IDCH3O2H,k) * xnew(OXYO2H)  &
    + rcphot(IDHCHO_H,k) *0.9 * xnew(C5DICAROOH)  &
    + rct(41,k) * xnew(TERPO2) * xnew(NO)  &
    + 2.3e-12 * xnew(TERPO2) * xnew(NO3)  &
    + 0.3* rct(41,k) * xnew(HMVKO2) * xnew(NO)  &
    + 0.3* 2.3e-12 * xnew(HMVKO2) * xnew(NO3)  &
    + 0.28* rct(99,k) * xnew(MVK) * xnew(O3)  &
    + rct(76,k) * xnew(APINENE) * xnew(O3)  &
    + AQRCK(ICLOHSO2,K)*2.00e-12 * xnew(OH) * xnew(SO2)
  
  L = rct(8,k) * xnew(O3)  &
    + rct(10,k) * xnew(NO)  &
    + rct(12,k) * xnew(OH)  &
    + rct(16,k) * xnew(HO2)  &
    + rct(16,k) * xnew(HO2)  &
    + rct(17,k) * xnew(HO2)  &
    + rct(17,k) * xnew(HO2)  &
    + 4.00e-12 * xnew(NO3)  &
    + rct(19,k) * xnew(NO2)  &
    + rct(27,k) * xnew(CH3O2)  &
    + rct(34,k) * xnew(C2H5O2)  &
    + rct(38,k) * xnew(CH3CO3)  &
    + rct(40,k) * xnew(C4H9O2)  &
    + rct(40,k) * xnew(MEKO2)  &
    + rct(42,k) * xnew(ETRO2)  &
    + rct(45,k) * xnew(PRRO2)  &
    + rct(49,k) * xnew(OXYO2)  &
    + rct(51,k) * xnew(C5DICARBO2)  &
    + rct(51,k) * xnew(ISRO2)  &
    + rct(40,k) * xnew(MACRO2)
  
  !-> HO2 cont.
    L = L +   rct(96,k) * xnew(TERPO2)  &
    + rct(40,k) * xnew(HMVKO2)  &
    + rct(104,k)
  
  xnew(HO2) = (xold(HO2) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3O2

  P = rct(22,k) * xnew(OH) * xnew(CH4)  &
    + 0.6* rct(29,k) * xnew(CH3OOH) * xnew(OH)  &
    + 0.05* rct(35,k) * xnew(OH) * xnew(CH3CHO)  &
    + rct(37,k) * xnew(CH3CO3) * xnew(NO)  &
    + 4.0e-12 * xnew(CH3CO3) * xnew(NO3)  &
    + 0.454* rct(38,k) * xnew(CH3CO3) * xnew(HO2)  &
    + 0.7* 1.0e-11 * xnew(RO2POOL) * xnew(CH3CO3)  &
    + 0.3* rct(44,k) * xnew(O3) * xnew(C3H6)  &
    + 0.08* rct(56,k) * xnew(O3) * xnew(C5H8)  &
    + rcphot(IDCH3CHO,k) * xnew(CH3CHO)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3CO3H)
  
  L = 1.2e-12 * xnew(NO3)  &
    + rct(24,k) * xnew(NO)  &
    + rct(25,k) * xnew(RO2POOL)  &
    + rct(26,k) * xnew(RO2POOL)  &
    + rct(27,k) * xnew(HO2)
  
  xnew(CH3O2) = (xold(CH3O2) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5O2

  P = rct(32,k) * xnew(OH) * xnew(C2H6)  &
    + 0.31* 1.17e-11 * xnew(C2H5OOH) * xnew(OH)  &
    + 0.264* rct(41,k) * xnew(NO) * xnew(C4H9O2)  &
    + 0.284* 2.3e-12 * xnew(NO3) * xnew(C4H9O2)  &
    + 0.348* rcphot(IDMEK,k) *0.78 * xnew(SC4H9NO3)  &
    + rct(46,k) * xnew(NO3) * xnew(C3H6)  &
    + rcphot(IDMEK,k) * xnew(MEK)  &
    + 0.284* rcphot(IDCH3O2H,k) * xnew(BURO2H)
  
  L = rct(33,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + rct(34,k) * xnew(HO2)
  
  xnew(C2H5O2) = (xold(C2H5O2) + dt2 * P) / (1.0 + dt2 * L)


!-> C4H9O2

  P = rct(39,k) * xnew(OH) * xnew(NC4H10)  &
    + 0.153* 2.4e-11 * xnew(BURO2H) * xnew(OH)
  
  L = rct(40,k) * xnew(HO2)  &
    + rct(41,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)
  
  xnew(C4H9O2) = (xold(C4H9O2) + dt2 * P) / (1.0 + dt2 * L)


!-> ISRO2

  P = rct(55,k) * xnew(OH) * xnew(C5H8)
  
  L = rct(41,k) * xnew(NO)  &
    + rct(51,k) * xnew(HO2)  &
    + 4.00e-12 * xnew(RO2POOL)  &
    + rct(65,k)  &
    + rct(66,k)
  
  xnew(ISRO2) = (xold(ISRO2) + dt2 * P) / (1.0 + dt2 * L)


!-> ETRO2

  P = 0.05* 3.22e-12 * xnew(OH) * xnew(C2H5OH)  &
    + 0.21* 1.7475e-11 * xnew(ETRO2H) * xnew(OH)  &
    + rct(72,k) * xnew(C2H4) * xnew(OH)
  
  L = rct(41,k) * xnew(NO)  &
    + rct(42,k) * xnew(HO2)
  
  xnew(ETRO2) = (xold(ETRO2) + dt2 * P) / (1.0 + dt2 * L)


!-> PRRO2

  P = 0.13* 2.81e-11 * xnew(PRRO2H) * xnew(OH)  &
    + rct(73,k) * xnew(OH) * xnew(C3H6)
  
  L = rct(41,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + rct(45,k) * xnew(HO2)
  
  xnew(PRRO2) = (xold(PRRO2) + dt2 * P) / (1.0 + dt2 * L)


!-> OXYO2

  P = 0.259* rct(47,k) * xnew(BENZENE) * xnew(OH)  &
    + 0.564* rct(48,k) * xnew(TOLUENE) * xnew(OH)  &
    + 0.55* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + 7.62e-11 * xnew(OXYO2H) * xnew(OH)
  
  L = rct(41,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + rct(49,k) * xnew(HO2)
  
  xnew(OXYO2) = (xold(OXYO2) + dt2 * P) / (1.0 + dt2 * L)


!-> MEKO2

  P = 1.1e-12 * xnew(OH) * xnew(MEK)  &
    + 0.174* 2.1e-11 * xnew(MEKO2H) * xnew(OH)
  
  L = rct(41,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + rct(40,k) * xnew(HO2)  &
    + 1.04e-12 * xnew(RO2POOL)
  
  xnew(MEKO2) = (xold(MEKO2) + dt2 * P) / (1.0 + dt2 * L)


!-> C5DICARBO2

  P = 6.2e-11 * xnew(C5DICARB) * xnew(OH)  &
    + rcphot(IDNO2,k) *0.2 * xnew(C5DICARB)  &
    + 0.03* 1.25e-10 * xnew(C5DICAROOH) * xnew(OH)
  
  L = rct(41,k) * xnew(NO)  &
    + rct(50,k) * xnew(NO3)  &
    + rct(51,k) * xnew(HO2)
  
  xnew(C5DICARBO2) = (xold(C5DICARBO2) + dt2 * P) / (1.0 + dt2 * L)


!-> MACRO2

  P = 0.1* rct(56,k) * xnew(O3) * xnew(C5H8)  &
    + 9.13e-12 * xnew(OH) * xnew(IEPOX)  &
    + rct(58,k) * xnew(OH) * xnew(MACR)  &
    + rct(59,k) * xnew(OH) * xnew(MACR)  &
    + rct(64,k) * xnew(MPAN)
  
  L = rct(41,k) * xnew(NO)  &
    + rct(40,k) * xnew(HO2)  &
    + rct(63,k) * xnew(NO2)
  
  xnew(MACRO2) = (xold(MACRO2) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3CO3

  P = 0.95* rct(35,k) * xnew(OH) * xnew(CH3CHO)  &
    + rct(36,k) * xnew(CH3CHO) * xnew(NO3)  &
    + 1.1e-11 * xnew(CH3CO3H) * xnew(OH)  &
    + rct(41,k) * xnew(MEKO2) * xnew(NO)  &
    + 2.3e-12 * xnew(MEKO2) * xnew(NO3)  &
    + 0.75* 1.04e-12 * xnew(RO2POOL) * xnew(MEKO2)  &
    + 2.3e-13 * xnew(BIACET) * xnew(OH)  &
    + 2* rcphot(IDCH3COY,k) * xnew(BIACET)  &
    + 1.72e-11 * xnew(C54CO) * xnew(OH)  &
    + rct(52,k) * xnew(C54CO) * xnew(NO3)  &
    + rcphot(IDRCOCHO,k) * xnew(C54CO)  &
    + rcphot(IDCH3COY,k) *2 * xnew(C54CO)  &
    + rct(54,k) * xnew(OH) * xnew(MGLYOX)  &
    + 0.1* rct(56,k) * xnew(O3) * xnew(C5H8)  &
    + 0.1* rct(60,k) * xnew(O3) * xnew(MACR)  &
    + 0.1* rct(61,k) * xnew(O3) * xnew(MACR)  &
    + 0.25* rct(41,k) * xnew(NO) * xnew(MACRO2)  &
    + rcphot(IDHCHO_H2,k) *0.0474 * xnew(MACR)  &
    + rcphot(IDMEK,k) * xnew(ACETOL)  &
    + rct(64,k) * xnew(PAN)
  
  !-> CH3CO3 cont.
    P = P +   rcphot(IDMEK,k) * xnew(MEK)  &
    + rcphot(IDRCOCHO,k) * xnew(MGLYOX)  &
    + rcphot(IDCH3O2H,k) * xnew(MEKO2H)  &
    + rcphot(IDMEK,k) * xnew(MEKO2H)  &
    + 0.7* rct(41,k) * xnew(HMVKO2) * xnew(NO)  &
    + 0.7* 2.3e-12 * xnew(HMVKO2) * xnew(NO3)  &
    + 0.28* rct(99,k) * xnew(MVK) * xnew(O3)
  
  L = rct(37,k) * xnew(NO)  &
    + 4.0e-12 * xnew(NO3)  &
    + rct(38,k) * xnew(HO2)  &
    + 1.0e-11 * xnew(RO2POOL)  &
    + rct(71,k) * xnew(NO2)
  
  xnew(CH3CO3) = (xold(CH3CO3) + dt2 * P) / (1.0 + dt2 * L)


!-> TERPO2

  P = rct(95,k) * xnew(APINENE) * xnew(OH)  &
    + rct(97,k) * xnew(TERPOOH) * xnew(OH)  &
    + rct(77,k) * xnew(APINENE) * xnew(NO3)
  
  L = rct(41,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + rct(96,k) * xnew(HO2)
  
  xnew(TERPO2) = (xold(TERPO2) + dt2 * P) / (1.0 + dt2 * L)


!-> HMVKO2

  P = rct(98,k) * xnew(MVK) * xnew(OH)
  
  L = rct(41,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + rct(40,k) * xnew(HO2)
  
  xnew(HMVKO2) = (xold(HMVKO2) + dt2 * P) / (1.0 + dt2 * L)


!-> RO2POOL

  ! P = 0.0
  
  ! L = 0.0
  
  xnew(RO2POOL) = sum(xnew(RO2_GROUP))


!-> NMVOC

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for NMVOC! xnew(NMVOC) = max(0.0, xold(NMVOC))


!-> O3

  P = rct(1,k) * xnew(OP)  &
    + 0.162* rct(38,k) * xnew(CH3CO3) * xnew(HO2)
  
  L = rct(5,k) * xnew(NO)  &
    + rct(6,k) * xnew(NO2)  &
    + rct(7,k) * xnew(OH)  &
    + rct(8,k) * xnew(HO2)  &
    + rct(43,k) * xnew(C2H4)  &
    + rct(44,k) * xnew(C3H6)  &
    + rct(56,k) * xnew(C5H8)  &
    + 9.20e-17 * xnew(ISON)  &
    + rct(60,k) * xnew(MACR)  &
    + rct(61,k) * xnew(MACR)  &
    + rcphot(IDO3_O1D,k)  &
    + rcphot(IDO3_O3P,k)  &
    + rct(99,k) * xnew(MVK)  &
    + rct(76,k) * xnew(APINENE)  &
    + AQRCK(ICLRC2,K) * xnew(SO2)
  
  xnew(O3) = (xold(O3) + dt2 * P) / (1.0 + dt2 * L)


!-> NO

  P = rct(11,k) * xnew(NO2) * xnew(NO3)  &
    + rcphot(IDNO2,k) * xnew(NO2)  &
    + 0.127* rcphot(IDNO3,k) * xnew(NO3)  &
    + rcphot(IDHONO,k) * xnew(HONO)  &
    + rcemis(NO,k)
  
  L = rct(5,k) * xnew(O3)  &
    + rct(9,k) * xnew(NO3)  &
    + rct(10,k) * xnew(HO2)  &
    + rct(24,k) * xnew(CH3O2)  &
    + rct(33,k) * xnew(C2H5O2)  &
    + rct(37,k) * xnew(CH3CO3)  &
    + rct(41,k) * xnew(C4H9O2)  &
    + rct(41,k) * xnew(MEKO2)  &
    + rct(41,k) * xnew(ETRO2)  &
    + rct(41,k) * xnew(PRRO2)  &
    + rct(41,k) * xnew(OXYO2)  &
    + rct(41,k) * xnew(C5DICARBO2)  &
    + rct(41,k) * xnew(ISRO2)  &
    + rct(41,k) * xnew(MACRO2)  &
    + rct(67,k) * xnew(OP)  &
    + rct(74,k) * xnew(OH)  &
    + rct(41,k) * xnew(TERPO2)  &
    + rct(41,k) * xnew(HMVKO2)
  
  xnew(NO) = (xold(NO) + dt2 * P) / (1.0 + dt2 * L)


!-> NO2

  P = rct(5,k) * xnew(O3) * xnew(NO)  &
    + rct(9,k) * xnew(NO) * xnew(NO3)  &
    + rct(9,k) * xnew(NO) * xnew(NO3)  &
    + rct(10,k) * xnew(NO) * xnew(HO2)  &
    + rct(11,k) * xnew(NO2) * xnew(NO3)  &
    + rct(18,k) * xnew(OH) * xnew(HONO)  &
    + 2.00e-11 * xnew(OH) * xnew(NO3)  &
    + 4.00e-12 * xnew(HO2) * xnew(NO3)  &
    + rct(20,k) * xnew(HO2NO2)  &
    + rct(21,k) * xnew(OH) * xnew(HO2NO2)  &
    + 1.2e-12 * xnew(CH3O2) * xnew(NO3)  &
    + rct(24,k) * xnew(CH3O2) * xnew(NO)  &
    + rct(33,k) * xnew(C2H5O2) * xnew(NO)  &
    + 2.3e-12 * xnew(C2H5O2) * xnew(NO3)  &
    + rct(37,k) * xnew(CH3CO3) * xnew(NO)  &
    + 4.0e-12 * xnew(CH3CO3) * xnew(NO3)  &
    + 0.917* rct(41,k) * xnew(NO) * xnew(C4H9O2)  &
    + 2.3e-12 * xnew(NO3) * xnew(C4H9O2)  &
    + 8.6e-13 * xnew(SC4H9NO3) * xnew(OH)  &
    + rcphot(IDMEK,k) *0.78 * xnew(SC4H9NO3)
  
  !-> NO2 cont.
    P = P +   rct(41,k) * xnew(MEKO2) * xnew(NO)  &
    + 2.3e-12 * xnew(MEKO2) * xnew(NO3)  &
    + rct(41,k) * xnew(ETRO2) * xnew(NO)  &
    + rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 2.3e-12 * xnew(NO3) * xnew(PRRO2)  &
    + rct(46,k) * xnew(NO3) * xnew(C3H6)  &
    + rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + 2.3e-12 * xnew(OXYO2) * xnew(NO3)  &
    + rct(41,k) * xnew(C5DICARBO2) * xnew(NO)  &
    + rct(50,k) * xnew(C5DICARBO2) * xnew(NO3)  &
    + 0.9* rct(41,k) * xnew(NO) * xnew(ISRO2)  &
    + 0.22* 3.34e-11 * xnew(OH) * xnew(ISON)  &
    + 0.33* 9.20e-17 * xnew(O3) * xnew(ISON)  &
    + rct(41,k) * xnew(NO) * xnew(MACRO2)  &
    + 2.90e-11 * xnew(OH) * xnew(MPAN)  &
    + 1.55e-12 * xnew(OH) * xnew(NALD)  &
    + rcphot(IDMEK,k) *0.78 * xnew(ISON)  &
    + rcphot(IDCH3CHO,k) * xnew(NALD)  &
    + rct(64,k) * xnew(MPAN)  &
    + rct(67,k) * xnew(OP) * xnew(NO)
  
  !-> NO2 cont.
    P = P +   rct(69,k) * xnew(N2O5)  &
    + rct(64,k) * xnew(PAN)  &
    + 0.873* rcphot(IDNO3,k) * xnew(NO3)  &
    + rcphot(IDHNO3,k) * xnew(HNO3)  &
    + 0.667* rcphot(IDHO2NO2,k) * xnew(HO2NO2)  &
    + rcemis(NO2,k)  &
    + rct(41,k) * xnew(TERPO2) * xnew(NO)  &
    + 2.3e-12 * xnew(TERPO2) * xnew(NO3)  &
    + rct(41,k) * xnew(HMVKO2) * xnew(NO)  &
    + 2.3e-12 * xnew(HMVKO2) * xnew(NO3)  &
    + rct(77,k) * xnew(APINENE) * xnew(NO3)
  
  L = rct(6,k) * xnew(O3)  &
    + rct(11,k) * xnew(NO3)  &
    + rct(19,k) * xnew(HO2)  &
    + rct(63,k) * xnew(MACRO2)  &
    + rct(68,k) * xnew(NO3)  &
    + rct(70,k) * xnew(OH)  &
    + rct(71,k) * xnew(CH3CO3)  &
    + rcphot(IDNO2,k)
  
  xnew(NO2) = (xold(NO2) + dt2 * P) / (1.0 + dt2 * L)


!-> NO3

  P = rct(6,k) * xnew(O3) * xnew(NO2)  &
    + rct(15,k) * xnew(OH) * xnew(HNO3)  &
    + rct(69,k) * xnew(N2O5)  &
    + 0.333* rcphot(IDHO2NO2,k) * xnew(HO2NO2)
  
  L = rct(9,k) * xnew(NO)  &
    + rct(11,k) * xnew(NO2)  &
    + 2.00e-11 * xnew(OH)  &
    + 4.00e-12 * xnew(HO2)  &
    + 1.2e-12 * xnew(CH3O2)  &
    + rct(31,k) * xnew(HCHO)  &
    + 2.3e-12 * xnew(C2H5O2)  &
    + rct(36,k) * xnew(CH3CHO)  &
    + 4.0e-12 * xnew(CH3CO3)  &
    + 2.3e-12 * xnew(C4H9O2)  &
    + 2.3e-12 * xnew(MEKO2)  &
    + 2.3e-12 * xnew(PRRO2)  &
    + rct(46,k) * xnew(C3H6)  &
    + 2.3e-12 * xnew(OXYO2)  &
    + rct(50,k) * xnew(C5DICARBO2)  &
    + rct(52,k) * xnew(C54CO)  &
    + rct(57,k) * xnew(C5H8)  &
    + rct(68,k) * xnew(NO2)  &
    + rcphot(IDNO3,k)  &
    + 2.3e-12 * xnew(TERPO2)
  
  !-> NO3 cont.
    L = L +   2.3e-12 * xnew(HMVKO2)  &
    + rct(77,k) * xnew(APINENE)  &
    + rct(101,k)
  
  xnew(NO3) = (xold(NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> N2O5

  P = rct(68,k) * xnew(NO2) * xnew(NO3)
  
  L = rct(69,k)  &
    + rct(100,k)
  
  xnew(N2O5) = (xold(N2O5) + dt2 * P) / (1.0 + dt2 * L)


!-> H2

  P = 0.18* rct(43,k) * xnew(C2H4) * xnew(O3)  &
    + 0.0684* rct(44,k) * xnew(O3) * xnew(C3H6)  &
    + rcphot(IDHCHO_H2,k) * xnew(HCHO)
  
  L = rct(14,k) * xnew(OH)
  
  xnew(H2) = (xold(H2) + dt2 * P) / (1.0 + dt2 * L)


!-> H2O2

  P = rct(16,k) * xnew(HO2) * xnew(HO2)  &
    + rct(17,k) * xnew(HO2) * xnew(HO2)  &
    + 0.1575* rct(43,k) * xnew(C2H4) * xnew(O3)  &
    + 0.169* rct(44,k) * xnew(O3) * xnew(C3H6)  &
    + 0.09* rct(56,k) * xnew(O3) * xnew(C5H8)  &
    + 0.09* rct(99,k) * xnew(MVK) * xnew(O3)  &
    + 0.5* rct(104,k) * xnew(HO2)
  
  L = rct(13,k) * xnew(OH)  &
    + rcphot(IDH2O2,k)  &
    + AQRCK(ICLRC1,K) * xnew(SO2)  &
    + 1.0e-5
  
  xnew(H2O2) = (xold(H2O2) + dt2 * P) / (1.0 + dt2 * L)


!-> HONO

  P = rct(74,k) * xnew(OH) * xnew(NO)
  
  L = rct(18,k) * xnew(OH)  &
    + rcphot(IDHONO,k)
  
  xnew(HONO) = (xold(HONO) + dt2 * P) / (1.0 + dt2 * L)


!-> HNO3

  P = rct(31,k) * xnew(NO3) * xnew(HCHO)  &
    + rct(36,k) * xnew(CH3CHO) * xnew(NO3)  &
    + rct(52,k) * xnew(C54CO) * xnew(NO3)  &
    + rct(70,k) * xnew(NO2) * xnew(OH)  &
    + 2* rct(100,k) * xnew(N2O5)  &
    + rct(101,k) * xnew(NO3)  &
    + rct(70,k) * xnew(shipNOx) * xnew(OH)  &
    + 3.2e-5 * xnew(shipNOx)
  
  L = rct(15,k) * xnew(OH)  &
    + rcphot(IDHNO3,k)  &
    + rct(102,k)  &
    + rct(103,k)
  
  xnew(HNO3) = (xold(HNO3) + dt2 * P) / (1.0 + dt2 * L)


!-> HO2NO2

  P = rct(19,k) * xnew(HO2) * xnew(NO2)
  
  L = rct(20,k)  &
    + rct(21,k) * xnew(OH)  &
    + rcphot(IDHO2NO2,k)
  
  xnew(HO2NO2) = (xold(HO2NO2) + dt2 * P) / (1.0 + dt2 * L)


!-> CO

  P = rct(30,k) * xnew(OH) * xnew(HCHO)  &
    + rct(31,k) * xnew(NO3) * xnew(HCHO)  &
    + 0.05* rct(35,k) * xnew(OH) * xnew(CH3CHO)  &
    + 2.3e-13 * xnew(BIACET) * xnew(OH)  &
    + 0.35* rct(43,k) * xnew(C2H4) * xnew(O3)  &
    + 0.483* rct(44,k) * xnew(O3) * xnew(C3H6)  &
    + 3* 1.72e-11 * xnew(C54CO) * xnew(OH)  &
    + 3* rct(52,k) * xnew(C54CO) * xnew(NO3)  &
    + 3* rcphot(IDRCOCHO,k) * xnew(C54CO)  &
    + 3* rcphot(IDCH3COY,k) *2 * xnew(C54CO)  &
    + 1.416* rct(53,k) * xnew(OH) * xnew(GLYOX)  &
    + rct(54,k) * xnew(OH) * xnew(MGLYOX)  &
    + 0.14* rct(56,k) * xnew(O3) * xnew(C5H8)  &
    + 0.22* rct(60,k) * xnew(O3) * xnew(MACR)  &
    + 0.22* rct(61,k) * xnew(O3) * xnew(MACR)  &
    + 0.25* rct(41,k) * xnew(NO) * xnew(MACRO2)  &
    + 3.00e-11 * xnew(OH) * xnew(MACROOH)  &
    + 1.55e-12 * xnew(OH) * xnew(NALD)  &
    + rcphot(IDHCHO_H2,k) *0.0474 * xnew(MACR)  &
    + 0.5* rcphot(IDCH3O2H,k) * xnew(MACROOH)
  
  !-> CO cont.
    P = P +   rcphot(IDCH3CHO,k) * xnew(NALD)  &
    + 7.61e-11 * xnew(OH) * xnew(HPALD)  &
    + rcphot(IDHCHO_H2,k) *24. * xnew(PACALD)  &
    + rcphot(IDHCHO_H,k) * xnew(HCHO)  &
    + rcphot(IDHCHO_H2,k) * xnew(HCHO)  &
    + rcphot(IDCH3CHO,k) * xnew(CH3CHO)  &
    + 1.9* rcphot(IDCHOCHO,k) * xnew(GLYOX)  &
    + rcphot(IDRCOCHO,k) * xnew(MGLYOX)  &
    + rcemis(CO,k)  &
    + rct(40,k) * xnew(HMVKO2) * xnew(HO2)  &
    + 0.56* rct(99,k) * xnew(MVK) * xnew(O3)
  
  L = rct(23,k) * xnew(OH)
  
  xnew(CO) = (xold(CO) + dt2 * P) / (1.0 + dt2 * L)


!-> CH4

  P = 0.04* rct(44,k) * xnew(O3) * xnew(C3H6)
  
  L = rct(22,k) * xnew(OH)
  
  xnew(CH4) = (xold(CH4) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H6

  P = 0.073* 1.04e-12 * xnew(RO2POOL) * xnew(MEKO2)  &
    + 0.0603* 2.1e-11 * xnew(MEKO2H) * xnew(OH)  &
    + rcemis(C2H6,k)  &
    + 0.19* rct(41,k) * xnew(TERPO2) * xnew(NO)  &
    + 0.19* 2.3e-12 * xnew(TERPO2) * xnew(NO3)
  
  L = rct(32,k) * xnew(OH)
  
  xnew(C2H6) = (xold(C2H6) + dt2 * P) / (1.0 + dt2 * L)


!-> NC4H10

  P = rcemis(NC4H10,k)
  
  L = rct(39,k) * xnew(OH)
  
  xnew(NC4H10) = (xold(NC4H10) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H4

  P = rcemis(C2H4,k)
  
  L = rct(43,k) * xnew(O3)  &
    + rct(72,k) * xnew(OH)
  
  xnew(C2H4) = (xold(C2H4) + dt2 * P) / (1.0 + dt2 * L)


!-> C3H6

  P = rcemis(C3H6,k)
  
  L = rct(44,k) * xnew(O3)  &
    + rct(46,k) * xnew(NO3)  &
    + rct(73,k) * xnew(OH)
  
  xnew(C3H6) = (xold(C3H6) + dt2 * P) / (1.0 + dt2 * L)


!-> BENZENE

  P = rcemis(BENZENE,k)
  
  L = rct(47,k) * xnew(OH)
  
  xnew(BENZENE) = (xold(BENZENE) + dt2 * P) / (1.0 + dt2 * L)


!-> TOLUENE

  P = rcemis(TOLUENE,k)
  
  L = rct(48,k) * xnew(OH)
  
  xnew(TOLUENE) = (xold(TOLUENE) + dt2 * P) / (1.0 + dt2 * L)


!-> OXYL

  P = rcemis(OXYL,k)
  
  L = 1.36e-11 * xnew(OH)
  
  xnew(OXYL) = (xold(OXYL) + dt2 * P) / (1.0 + dt2 * L)


!-> C5H8

  P = rcemis(C5H8,k)  &
    + RCBIO(1,K)
  
  L = rct(55,k) * xnew(OH)  &
    + rct(56,k) * xnew(O3)  &
    + rct(57,k) * xnew(NO3)
  
  xnew(C5H8) = (xold(C5H8) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3OH

  P = 0.5* rct(26,k) * xnew(RO2POOL) * xnew(CH3O2)  &
    + 0.05* rct(44,k) * xnew(O3) * xnew(C3H6)  &
    + rcemis(CH3OH,k)
  
  L = rct(28,k) * xnew(OH)
  
  xnew(CH3OH) = (xold(CH3OH) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5OH

  P = rcemis(C2H5OH,k)
  
  L = 3.22e-12 * xnew(OH)
  
  xnew(C2H5OH) = (xold(C2H5OH) + dt2 * P) / (1.0 + dt2 * L)


!-> HCHO

  P = 1.2e-12 * xnew(CH3O2) * xnew(NO3)  &
    + rct(24,k) * xnew(CH3O2) * xnew(NO)  &
    + rct(25,k) * xnew(RO2POOL) * xnew(CH3O2)  &
    + 0.5* rct(26,k) * xnew(RO2POOL) * xnew(CH3O2)  &
    + 0.1* rct(27,k) * xnew(CH3O2) * xnew(HO2)  &
    + rct(28,k) * xnew(CH3OH) * xnew(OH)  &
    + 0.4* rct(29,k) * xnew(CH3OOH) * xnew(OH)  &
    + 0.1095* 1.04e-12 * xnew(RO2POOL) * xnew(MEKO2)  &
    + 2.3e-13 * xnew(BIACET) * xnew(OH)  &
    + 2* rct(41,k) * xnew(ETRO2) * xnew(NO)  &
    + 1.1575* rct(43,k) * xnew(C2H4) * xnew(O3)  &
    + 0.68* rct(44,k) * xnew(O3) * xnew(C3H6)  &
    + rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 2.3e-12 * xnew(NO3) * xnew(PRRO2)  &
    + rct(46,k) * xnew(NO3) * xnew(C3H6)  &
    + 0.58* rct(56,k) * xnew(O3) * xnew(C5H8)  &
    + 0.9* rct(41,k) * xnew(NO) * xnew(ISRO2)  &
    + 4.00e-12 * xnew(RO2POOL) * xnew(ISRO2)  &
    + 0.22* 3.34e-11 * xnew(OH) * xnew(ISON)  &
    + 0.33* 9.20e-17 * xnew(O3) * xnew(ISON)
  
  !-> HCHO cont.
    P = P +   0.75* rct(41,k) * xnew(NO) * xnew(MACRO2)  &
    + 1.55e-12 * xnew(OH) * xnew(NALD)  &
    + rcphot(IDCH3O2H,k) * xnew(ISRO2H)  &
    + rcphot(IDMEK,k) *0.78 * xnew(ISON)  &
    + rcphot(IDHCHO_H2,k) *0.0474 * xnew(MACR)  &
    + 0.5* rcphot(IDCH3O2H,k) * xnew(MACROOH)  &
    + rcphot(IDMEK,k) * xnew(ACETOL)  &
    + rcphot(IDCH3CHO,k) * xnew(NALD)  &
    + rct(66,k) * xnew(ISRO2)  &
    + 7.61e-11 * xnew(OH) * xnew(HPALD)  &
    + 0.1* rcphot(IDCHOCHO,k) * xnew(GLYOX)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3OOH)  &
    + 1.56* rcphot(IDCH3O2H,k) * xnew(ETRO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(PRRO2H)  &
    + rcemis(HCHO,k)  &
    + 0.3* rct(41,k) * xnew(HMVKO2) * xnew(NO)  &
    + 0.3* 2.3e-12 * xnew(HMVKO2) * xnew(NO3)  &
    + 0.645* rct(99,k) * xnew(MVK) * xnew(O3)
  
  L = rct(30,k) * xnew(OH)  &
    + rct(31,k) * xnew(NO3)  &
    + rcphot(IDHCHO_H,k)  &
    + rcphot(IDHCHO_H2,k)
  
  xnew(HCHO) = (xold(HCHO) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3CHO

  P = rct(33,k) * xnew(C2H5O2) * xnew(NO)  &
    + 2.3e-12 * xnew(C2H5O2) * xnew(NO3)  &
    + 0.69* 1.17e-11 * xnew(C2H5OOH) * xnew(OH)  &
    + 0.95* 3.22e-12 * xnew(OH) * xnew(C2H5OH)  &
    + 0.264* rct(41,k) * xnew(NO) * xnew(C4H9O2)  &
    + 0.284* 2.3e-12 * xnew(NO3) * xnew(C4H9O2)  &
    + 0.348* rcphot(IDMEK,k) *0.78 * xnew(SC4H9NO3)  &
    + rct(41,k) * xnew(MEKO2) * xnew(NO)  &
    + 2.3e-12 * xnew(MEKO2) * xnew(NO3)  &
    + 0.6405* 1.04e-12 * xnew(RO2POOL) * xnew(MEKO2)  &
    + 0.79* 1.7475e-11 * xnew(ETRO2H) * xnew(OH)  &
    + 0.489* rct(44,k) * xnew(O3) * xnew(C3H6)  &
    + rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 2.3e-12 * xnew(NO3) * xnew(PRRO2)  &
    + rcphot(IDCH3O2H,k) * xnew(C2H5OOH)  &
    + 0.22* rcphot(IDCH3O2H,k) * xnew(ETRO2H)  &
    + 0.284* rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(PRRO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(MEKO2H)  &
    + rcphot(IDMEK,k) * xnew(MEKO2H)
  
  !-> CH3CHO cont.
    P = P +   rcemis(CH3CHO,k)  &
    + 0.7* rct(41,k) * xnew(HMVKO2) * xnew(NO)  &
    + 0.7* 2.3e-12 * xnew(HMVKO2) * xnew(NO3)  &
    + 0.1* rct(99,k) * xnew(MVK) * xnew(O3)
  
  L = rct(35,k) * xnew(OH)  &
    + rct(36,k) * xnew(NO3)  &
    + rcphot(IDCH3CHO,k)
  
  xnew(CH3CHO) = (xold(CH3CHO) + dt2 * P) / (1.0 + dt2 * L)


!-> MACR

  P = 0.65* rct(56,k) * xnew(O3) * xnew(C5H8)  &
    + 0.9* rct(41,k) * xnew(NO) * xnew(ISRO2)  &
    + 4.00e-12 * xnew(RO2POOL) * xnew(ISRO2)  &
    + 0.1* 8.94e-11 * xnew(OH) * xnew(ISRO2H)  &
    + 0.22* 3.34e-11 * xnew(OH) * xnew(ISON)  &
    + 0.33* 9.20e-17 * xnew(O3) * xnew(ISON)  &
    + rcphot(IDCH3O2H,k) * xnew(ISRO2H)  &
    + rcphot(IDMEK,k) *0.78 * xnew(ISON)  &
    + rct(66,k) * xnew(ISRO2)  &
    + rct(41,k) * xnew(TERPO2) * xnew(NO)  &
    + 2.3e-12 * xnew(TERPO2) * xnew(NO3)  &
    + rct(76,k) * xnew(APINENE) * xnew(O3)
  
  L = rct(58,k) * xnew(OH)  &
    + rct(59,k) * xnew(OH)  &
    + rct(60,k) * xnew(O3)  &
    + rct(61,k) * xnew(O3)  &
    + rcphot(IDHCHO_H2,k) *0.0474
  
  xnew(MACR) = (xold(MACR) + dt2 * P) / (1.0 + dt2 * L)


!-> MEK

  P = 0.847* 2.4e-11 * xnew(BURO2H) * xnew(OH)  &
    + 0.654* rct(41,k) * xnew(NO) * xnew(C4H9O2)  &
    + 0.716* 2.3e-12 * xnew(NO3) * xnew(C4H9O2)  &
    + 8.6e-13 * xnew(SC4H9NO3) * xnew(OH)  &
    + 0.652* rcphot(IDMEK,k) *0.78 * xnew(SC4H9NO3)  &
    + 0.716* rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + rcemis(MEK,k)
  
  L = 1.1e-12 * xnew(OH)  &
    + rcphot(IDMEK,k)
  
  xnew(MEK) = (xold(MEK) + dt2 * P) / (1.0 + dt2 * L)


!-> ACETOL

  P = 0.87* 2.81e-11 * xnew(PRRO2H) * xnew(OH)  &
    + 0.78* 3.34e-11 * xnew(OH) * xnew(ISON)  &
    + 0.25* rct(41,k) * xnew(NO) * xnew(MACRO2)  &
    + 2.90e-11 * xnew(OH) * xnew(MPAN)  &
    + 3.00e-11 * xnew(OH) * xnew(MACROOH)  &
    + 0.5* rcphot(IDCH3O2H,k) * xnew(MACROOH)  &
    + rct(40,k) * xnew(HMVKO2) * xnew(HO2)
  
  L = rct(62,k) * xnew(OH)  &
    + rcphot(IDMEK,k)
  
  xnew(ACETOL) = (xold(ACETOL) + dt2 * P) / (1.0 + dt2 * L)


!-> GLYOX

  P = rct(41,k) * xnew(C5DICARBO2) * xnew(NO)  &
    + rct(50,k) * xnew(C5DICARBO2) * xnew(NO3)  &
    + rcphot(IDHCHO_H,k) *0.9 * xnew(C5DICAROOH)  &
    + rcemis(GLYOX,k)
  
  L = rct(53,k) * xnew(OH)  &
    + rcphot(IDCHOCHO,k)
  
  xnew(GLYOX) = (xold(GLYOX) + dt2 * P) / (1.0 + dt2 * L)


!-> MGLYOX

  P = 0.0365* 1.04e-12 * xnew(RO2POOL) * xnew(MEKO2)  &
    + 0.1206* 2.1e-11 * xnew(MEKO2H) * xnew(OH)  &
    + rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + 2.3e-12 * xnew(OXYO2) * xnew(NO3)  &
    + rct(41,k) * xnew(C5DICARBO2) * xnew(NO)  &
    + rct(50,k) * xnew(C5DICARBO2) * xnew(NO3)  &
    + 0.9* rct(60,k) * xnew(O3) * xnew(MACR)  &
    + 0.9* rct(61,k) * xnew(O3) * xnew(MACR)  &
    + 0.5* rct(41,k) * xnew(NO) * xnew(MACRO2)  &
    + rct(62,k) * xnew(OH) * xnew(ACETOL)  &
    + 0.5* rcphot(IDCH3O2H,k) * xnew(MACROOH)  &
    + 7.61e-11 * xnew(OH) * xnew(HPALD)  &
    + rcphot(IDHCHO_H2,k) *24. * xnew(PACALD)  &
    + rcphot(IDCH3O2H,k) * xnew(OXYO2H)  &
    + rcphot(IDHCHO_H,k) *0.9 * xnew(C5DICAROOH)  &
    + rcemis(MGLYOX,k)  &
    + 0.3* rct(41,k) * xnew(HMVKO2) * xnew(NO)  &
    + 0.3* 2.3e-12 * xnew(HMVKO2) * xnew(NO3)  &
    + 0.62* rct(99,k) * xnew(MVK) * xnew(O3)
  
  L = rct(54,k) * xnew(OH)  &
    + rcphot(IDRCOCHO,k)
  
  xnew(MGLYOX) = (xold(MGLYOX) + dt2 * P) / (1.0 + dt2 * L)


!-> BIACET

  P = 0.2135* 1.04e-12 * xnew(RO2POOL) * xnew(MEKO2)  &
    + 0.7054* 2.1e-11 * xnew(MEKO2H) * xnew(OH)
  
  L = 2.3e-13 * xnew(OH)  &
    + rcphot(IDCH3COY,k)
  
  xnew(BIACET) = (xold(BIACET) + dt2 * P) / (1.0 + dt2 * L)


!-> C5DICARB

  P = rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + 2.3e-12 * xnew(OXYO2) * xnew(NO3)  &
    + rcphot(IDCH3O2H,k) * xnew(OXYO2H)
  
  L = 6.2e-11 * xnew(OH)  &
    + rcphot(IDNO2,k) *0.2
  
  xnew(C5DICARB) = (xold(C5DICARB) + dt2 * P) / (1.0 + dt2 * L)


!-> C5134CO2OH

  P = 0.97* 1.25e-10 * xnew(C5DICAROOH) * xnew(OH)
  
  L = 7.48e-11 * xnew(OH)
  
  xnew(C5134CO2OH) = (xold(C5134CO2OH) + dt2 * P) / (1.0 + dt2 * L)


!-> C54CO

  P = 7.48e-11 * xnew(C5134CO2OH) * xnew(OH)
  
  L = 1.72e-11 * xnew(OH)  &
    + rct(52,k) * xnew(NO3)  &
    + rcphot(IDRCOCHO,k)  &
    + rcphot(IDCH3COY,k) *2
  
  xnew(C54CO) = (xold(C54CO) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3OOH

  P = 0.9* rct(27,k) * xnew(CH3O2) * xnew(HO2)
  
  L = rct(29,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)  &
    + 1.0e-5
  
  xnew(CH3OOH) = (xold(CH3OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5OOH

  P = rct(34,k) * xnew(C2H5O2) * xnew(HO2)
  
  L = 1.17e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(C2H5OOH) = (xold(C2H5OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> BURO2H

  P = rct(40,k) * xnew(C4H9O2) * xnew(HO2)
  
  L = 2.4e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(BURO2H) = (xold(BURO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> ETRO2H

  P = rct(42,k) * xnew(ETRO2) * xnew(HO2)
  
  L = 1.7475e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(ETRO2H) = (xold(ETRO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> PRRO2H

  P = rct(45,k) * xnew(PRRO2) * xnew(HO2)
  
  L = 2.81e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(PRRO2H) = (xold(PRRO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> MEKO2H

  P = rct(40,k) * xnew(MEKO2) * xnew(HO2)
  
  L = 2.1e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)  &
    + rcphot(IDMEK,k)
  
  xnew(MEKO2H) = (xold(MEKO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> ISRO2H

  P = rct(51,k) * xnew(HO2) * xnew(ISRO2)
  
  L = 8.94e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(ISRO2H) = (xold(ISRO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> C5DICAROOH

  P = rct(51,k) * xnew(C5DICARBO2) * xnew(HO2)
  
  L = 1.25e-10 * xnew(OH)  &
    + rcphot(IDHCHO_H,k) *0.9
  
  xnew(C5DICAROOH) = (xold(C5DICAROOH) + dt2 * P) / (1.0 + dt2 * L)


!-> HPALD

  P = rct(65,k) * xnew(ISRO2)
  
  L = rcphot(IDHCHO_H2,k) *12.  &
    + 7.61e-11 * xnew(OH)
  
  xnew(HPALD) = (xold(HPALD) + dt2 * P) / (1.0 + dt2 * L)


!-> MACROOH

  P = rct(40,k) * xnew(HO2) * xnew(MACRO2)
  
  L = 3.00e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(MACROOH) = (xold(MACROOH) + dt2 * P) / (1.0 + dt2 * L)


!-> OXYO2H

  P = rct(49,k) * xnew(OXYO2) * xnew(HO2)
  
  L = 7.62e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(OXYO2H) = (xold(OXYO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3CO3H

  P = 0.384* rct(38,k) * xnew(CH3CO3) * xnew(HO2)
  
  L = 1.1e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(CH3CO3H) = (xold(CH3CO3H) + dt2 * P) / (1.0 + dt2 * L)


!-> PACALD

  P = rcphot(IDHCHO_H2,k) *12. * xnew(HPALD)
  
  L = rcphot(IDHCHO_H2,k) *24.
  
  xnew(PACALD) = (xold(PACALD) + dt2 * P) / (1.0 + dt2 * L)


!-> IEPOX

  P = 0.9* 8.94e-11 * xnew(OH) * xnew(ISRO2H)
  
  L = 9.13e-12 * xnew(OH)
  
  xnew(IEPOX) = (xold(IEPOX) + dt2 * P) / (1.0 + dt2 * L)


!-> SC4H9NO3

  P = 0.083* rct(41,k) * xnew(NO) * xnew(C4H9O2)
  
  L = 8.6e-13 * xnew(OH)  &
    + rcphot(IDMEK,k) *0.78
  
  xnew(SC4H9NO3) = (xold(SC4H9NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> NALD

  P = 0.78* 3.34e-11 * xnew(OH) * xnew(ISON)  &
    + 0.67* 9.20e-17 * xnew(O3) * xnew(ISON)
  
  L = 1.55e-12 * xnew(OH)  &
    + rcphot(IDCH3CHO,k)
  
  xnew(NALD) = (xold(NALD) + dt2 * P) / (1.0 + dt2 * L)


!-> ISON

  P = rct(57,k) * xnew(NO3) * xnew(C5H8)  &
    + 0.1* rct(41,k) * xnew(NO) * xnew(ISRO2)
  
  L = 3.34e-11 * xnew(OH)  &
    + 9.20e-17 * xnew(O3)  &
    + rcphot(IDMEK,k) *0.78
  
  xnew(ISON) = (xold(ISON) + dt2 * P) / (1.0 + dt2 * L)


!-> PAN

  P = rct(71,k) * xnew(CH3CO3) * xnew(NO2)
  
  L = rct(64,k)
  
  xnew(PAN) = (xold(PAN) + dt2 * P) / (1.0 + dt2 * L)


!-> MPAN

  P = rct(63,k) * xnew(MACRO2) * xnew(NO2)
  
  L = 2.90e-11 * xnew(OH)  &
    + rct(64,k)
  
  xnew(MPAN) = (xold(MPAN) + dt2 * P) / (1.0 + dt2 * L)


!-> SO2

  P = rcemis(SO2,k)
  
  L = AQRCK(ICLOHSO2,K)*2.00e-12 * xnew(OH)  &
    + AQRCK(ICLRC1,K) * xnew(H2O2)  &
    + AQRCK(ICLRC2,K) * xnew(O3)  &
    + AQRCK(ICLRC3,K)
  
  xnew(SO2) = (xold(SO2) + dt2 * P) / (1.0 + dt2 * L)


!-> TERPPeroxy

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for TERPPeroxy! xnew(TERPPeroxy) = max(0.0, xold(TERPPeroxy))


!-> VBS_TEST

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for VBS_TEST! xnew(VBS_TEST) = max(0.0, xold(VBS_TEST))


!-> APINENE

  P = RCBIO(2,K)
  
  L = rct(95,k) * xnew(OH)  &
    + rct(76,k) * xnew(O3)  &
    + rct(77,k) * xnew(NO3)
  
  xnew(APINENE) = (xold(APINENE) + dt2 * P) / (1.0 + dt2 * L)


!-> SQT_SOA_nv

  P = 0.00383* RCBIO(2,K)
  
  ! L = 0.0
  
  xnew(SQT_SOA_nv) = xold(SQT_SOA_nv) + dt2 * P


!-> TERPOOH

  P = rct(96,k) * xnew(TERPO2) * xnew(HO2)
  
  L = rct(97,k) * xnew(OH)
  
  xnew(TERPOOH) = (xold(TERPOOH) + dt2 * P) / (1.0 + dt2 * L)


!-> MVK

  P = rct(41,k) * xnew(TERPO2) * xnew(NO)  &
    + 2.3e-12 * xnew(TERPO2) * xnew(NO3)  &
    + rct(76,k) * xnew(APINENE) * xnew(O3)
  
  L = rct(98,k) * xnew(OH)  &
    + rct(99,k) * xnew(O3)
  
  xnew(MVK) = (xold(MVK) + dt2 * P) / (1.0 + dt2 * L)


!-> shipNOx

  P = rcemis(SHIPNOX,k)
  
  L = rct(70,k) * xnew(OH)  &
    + 3.2e-5
  
  xnew(shipNOx) = (xold(shipNOx) + dt2 * P) / (1.0 + dt2 * L)


!-> Dust_road_f

  P = rcemis(DUST_ROAD_F,K)
  
  ! L = 0.0
  
  xnew(Dust_road_f) = xold(Dust_road_f) + dt2 * P


!-> Dust_road_c

  P = rcemis(DUST_ROAD_C,K)
  
  ! L = 0.0
  
  xnew(Dust_road_c) = xold(Dust_road_c) + dt2 * P


!-> Dust_wb_f

  P = rcemis(DUST_WB_F,K)
  
  ! L = 0.0
  
  xnew(Dust_wb_f) = xold(Dust_wb_f) + dt2 * P


!-> Dust_wb_c

  P = rcemis(DUST_WB_C,K)
  
  ! L = 0.0
  
  xnew(Dust_wb_c) = xold(Dust_wb_c) + dt2 * P


!-> Dust_sah_f

  P = rcemis(DUST_SAH_F,K)
  
  ! L = 0.0
  
  xnew(Dust_sah_f) = xold(Dust_sah_f) + dt2 * P


!-> Dust_sah_c

  P = rcemis(DUST_SAH_C,K)
  
  ! L = 0.0
  
  xnew(Dust_sah_c) = xold(Dust_sah_c) + dt2 * P


