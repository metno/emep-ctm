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
    + rct(10,k) * xnew(NO)
  
  xnew(OP) = (xold(OP) + dt2 * P) / (1.0 + dt2 * L)


!-> OH

  P = 2* rct(4,k) * xnew(OD)  &
    + rct(14,k) * xnew(O3) * xnew(HO2)  &
    + rct(16,k) * xnew(NO) * xnew(HO2)  &
    + 4.0e-12 * xnew(HO2) * xnew(NO3)  &
    + 0.4* rct(34,k) * xnew(CH3OOH) * xnew(OH)  &
    + 0.69* 1.17e-11 * xnew(C2H5OOH) * xnew(OH)  &
    + 0.454* rct(45,k) * xnew(CH3CO3) * xnew(HO2)  &
    + 0.847* 2.4e-11 * xnew(BURO2H) * xnew(OH)  &
    + 0.826* 2.1e-11 * xnew(MEKO2H) * xnew(OH)  &
    + 0.79* 1.7475e-11 * xnew(ETRO2H) * xnew(OH)  &
    + 0.17* rct(51,k) * xnew(C2H4) * xnew(O3)  &
    + 0.3646* rct(53,k) * xnew(C3H6) * xnew(O3)  &
    + 0.584* rct(62,k) * xnew(GLYOX) * xnew(OH)  &
    + 0.25* rct(65,k) * xnew(C5H8) * xnew(O3)  &
    + 8.94e-11 * xnew(ISRO2H) * xnew(OH)  &
    + 0.67* 9.2e-17 * xnew(ISON) * xnew(O3)  &
    + 0.19* rct(69,k) * xnew(MACR) * xnew(O3)  &
    + 0.19* rct(70,k) * xnew(MACR) * xnew(O3)  &
    + 3.0e-11 * xnew(MACROOH) * xnew(OH)  &
    + rct(74,k) * xnew(ISRO2)
  
  !-> OH cont.
    P = P +   7.61e-11 * xnew(HPALD) * xnew(OH)  &
    + 0.7* rct(78,k) * xnew(APINENE) * xnew(O3)  &
    + 2* rcphot(IDH2O2,k) * xnew(H2O2)  &
    + rcphot(IDHONO,k) * xnew(HONO)  &
    + rcphot(IDHNO3,k) * xnew(HNO3)  &
    + 0.333* rcphot(IDHO2NO2,k) * xnew(HO2NO2)  &
    + rcphot(IDHCHO_H2,k) *12. * xnew(HPALD)  &
    + rcphot(IDHCHO_H2,k) *24. * xnew(PACALD)  &
    + rcphot(IDHCHO_H,k) *0.9 * xnew(C5DICAROOH)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3OOH)  &
    + rcphot(IDCH3O2H,k) *0.5 * xnew(C2H5OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(ETRO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3CO3H)  &
    + rcphot(IDCH3O2H,k) * xnew(MEKO2H)  &
    + rcphot(IDMEK,k) * xnew(MEKO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(MACROOH)  &
    + rcphot(IDCH3O2H,k) * xnew(OXYO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(ISRO2H)
  
  L = rct(5,k) * xnew(HO2)  &
    + rct(6,k) * xnew(H2O2)  &
    + rct(7,k)  &
    + rct(13,k) * xnew(O3)  &
    + rct(20,k) * xnew(NO2)  &
    + rct(21,k) * xnew(HNO3)  &
    + rct(22,k) * xnew(NO)  &
    + rct(23,k) * xnew(HONO)  &
    + 2.0e-11 * xnew(NO3)  &
    + rct(26,k) * xnew(HO2NO2)  &
    + rct(27,k)  &
    + rct(28,k) * xnew(CO)  &
    + rct(33,k) * xnew(CH3OH)  &
    + rct(34,k) * xnew(CH3OOH)  &
    + rct(35,k) * xnew(HCHO)  &
    + rct(37,k) * xnew(C2H6)  &
    + 1.17e-11 * xnew(C2H5OOH)  &
    + rct(40,k) * xnew(CH3CHO)  &
    + 1.1e-11 * xnew(CH3CO3H)  &
    + 3.22e-12 * xnew(C2H5OH)
  
  !-> OH cont.
    L = L +   rct(46,k) * xnew(NC4H10)  &
    + 2.4e-11 * xnew(BURO2H)  &
    + 8.6e-13 * xnew(SC4H9NO3)  &
    + 1.1e-12 * xnew(MEK)  &
    + 2.1e-11 * xnew(MEKO2H)  &
    + 2.3e-13 * xnew(BIACET)  &
    + rct(49,k) * xnew(C2H4)  &
    + 1.7475e-11 * xnew(ETRO2H)  &
    + rct(52,k) * xnew(C3H6)  &
    + rct(56,k) * xnew(BENZENE)  &
    + rct(57,k) * xnew(TOLUENE)  &
    + 1.36e-11 * xnew(OXYL)  &
    + 7.62e-11 * xnew(OXYO2H)  &
    + 6.2e-11 * xnew(C5DICARB)  &
    + rct(61,k) * xnew(C5DICAROOH)  &
    + rct(62,k) * xnew(GLYOX)  &
    + rct(63,k) * xnew(MGLYOX)  &
    + rct(64,k) * xnew(C5H8)  &
    + 8.94e-11 * xnew(ISRO2H)  &
    + 9.13e-12 * xnew(IEPOX)
  
  !-> OH cont.
    L = L +   3.34e-11 * xnew(ISON)  &
    + rct(67,k) * xnew(MACR)  &
    + rct(68,k) * xnew(MACR)  &
    + 2.9e-11 * xnew(MPAN)  &
    + 3.0e-11 * xnew(MACROOH)  &
    + rct(71,k) * xnew(ACETOL)  &
    + 1.55e-12 * xnew(NALD)  &
    + 7.61e-11 * xnew(HPALD)  &
    + rct(75,k) * xnew(APINENE)  &
    + rct(77,k) * xnew(TERPOOH)  &
    + AQRCK(ICLOHSO2,K)*2.00e-12 * xnew(SO2)  &
    + rct(20,k) * xnew(shipNOx)
  
  xnew(OH) = (xold(OH) + dt2 * P) / (1.0 + dt2 * L)


!-> HO2

  P = rct(6,k) * xnew(OH) * xnew(H2O2)  &
    + rct(7,k) * xnew(OH)  &
    + rct(13,k) * xnew(O3) * xnew(OH)  &
    + 2.0e-11 * xnew(OH) * xnew(NO3)  &
    + rct(25,k) * xnew(HO2NO2)  &
    + rct(28,k) * xnew(CO) * xnew(OH)  &
    + 1.2e-12 * xnew(CH3O2) * xnew(NO3)  &
    + rct(29,k) * xnew(CH3O2) * xnew(NO)  &
    + rct(30,k) * xnew(RO2POOL) * xnew(CH3O2)  &
    + rct(33,k) * xnew(CH3OH) * xnew(OH)  &
    + rct(35,k) * xnew(HCHO) * xnew(OH)  &
    + rct(36,k) * xnew(HCHO) * xnew(NO3)  &
    + rct(38,k) * xnew(C2H5O2) * xnew(NO)  &
    + 2.3e-12 * xnew(C2H5O2) * xnew(NO3)  &
    + 0.95* 3.22e-12 * xnew(C2H5OH) * xnew(OH)  &
    + 0.654* rct(48,k) * xnew(C4H9O2) * xnew(NO)  &
    + 0.716* 2.3e-12 * xnew(C4H9O2) * xnew(NO3)  &
    + rct(48,k) * xnew(ETRO2) * xnew(NO)  &
    + 0.27* rct(51,k) * xnew(C2H4) * xnew(O3)  &
    + 0.038* rct(53,k) * xnew(C3H6) * xnew(O3)
  
  !-> HO2 cont.
    P = P +   rct(48,k) * xnew(PRRO2) * xnew(NO)  &
    + 2.3e-12 * xnew(PRRO2) * xnew(NO3)  &
    + rct(48,k) * xnew(OXYO2) * xnew(NO)  &
    + 2.3e-12 * xnew(OXYO2) * xnew(NO3)  &
    + rct(48,k) * xnew(C5DICARBO2) * xnew(NO)  &
    + rct(59,k) * xnew(C5DICARBO2) * xnew(NO3)  &
    + 0.416* rct(62,k) * xnew(GLYOX) * xnew(OH)  &
    + 0.25* rct(65,k) * xnew(C5H8) * xnew(O3)  &
    + rct(66,k) * xnew(C5H8) * xnew(NO3)  &
    + 0.9* rct(48,k) * xnew(ISRO2) * xnew(NO)  &
    + 4.0e-12 * xnew(RO2POOL) * xnew(ISRO2)  &
    + 0.78* 3.34e-11 * xnew(ISON) * xnew(OH)  &
    + 0.32* rct(69,k) * xnew(MACR) * xnew(O3)  &
    + 0.32* rct(70,k) * xnew(MACR) * xnew(O3)  &
    + 0.75* rct(48,k) * xnew(MACRO2) * xnew(NO)  &
    + rct(71,k) * xnew(ACETOL) * xnew(OH)  &
    + rct(73,k) * xnew(ISRO2)  &
    + rct(48,k) * xnew(TERPO2) * xnew(NO)  &
    + 2.3e-12 * xnew(TERPO2) * xnew(NO3)  &
    + rct(78,k) * xnew(APINENE) * xnew(O3)
  
  !-> HO2 cont.
    P = P +   0.667* rcphot(IDHO2NO2,k) * xnew(HO2NO2)  &
    + 2* rcphot(IDHCHO_H,k) * xnew(HCHO)  &
    + rcphot(IDCH3CHO,k) * xnew(CH3CHO)  &
    + rcphot(IDCH3CHO,k) * xnew(NALD)  &
    + rcphot(IDMEK,k) * xnew(ACETOL)  &
    + rcphot(IDRCOCHO,k) * xnew(MGLYOX)  &
    + rcphot(IDHCHO_H2,k) *0.0474 * xnew(MACR)  &
    + rcphot(IDNO2,k) *0.2 * xnew(C5DICARB)  &
    + rcphot(IDHCHO_H2,k) *12. * xnew(HPALD)  &
    + rcphot(IDHCHO_H2,k) *24. * xnew(PACALD)  &
    + rcphot(IDHCHO_H,k) *0.9 * xnew(C5DICAROOH)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3OOH)  &
    + rcphot(IDCH3O2H,k) *0.5 * xnew(C2H5OOH)  &
    + 0.716* rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(ETRO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(MACROOH)  &
    + rcphot(IDCH3O2H,k) * xnew(OXYO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(ISRO2H)  &
    + 0.652* rcphot(IDMEK,k) *0.78 * xnew(SC4H9NO3)  &
    + rcphot(IDMEK,k) *0.78 * xnew(ISON)
  
  !-> HO2 cont.
    P = P +   2* rcphot(IDGLYOXA,k) * xnew(GLYOX)  &
    + AQRCK(ICLOHSO2,K)*2.00e-12 * xnew(OH) * xnew(SO2)
  
  L = rct(5,k) * xnew(OH)  &
    + rct(8,k) * xnew(HO2)  &
    + rct(8,k) * xnew(HO2)  &
    + rct(9,k) * xnew(HO2)  &
    + rct(9,k) * xnew(HO2)  &
    + rct(14,k) * xnew(O3)  &
    + rct(16,k) * xnew(NO)  &
    + 4.0e-12 * xnew(NO3)  &
    + rct(24,k) * xnew(NO2)  &
    + rct(32,k) * xnew(CH3O2)  &
    + rct(39,k) * xnew(C2H5O2)  &
    + rct(45,k) * xnew(CH3CO3)  &
    + rct(47,k) * xnew(C4H9O2)  &
    + rct(47,k) * xnew(MEKO2)  &
    + rct(50,k) * xnew(ETRO2)  &
    + rct(54,k) * xnew(PRRO2)  &
    + rct(58,k) * xnew(OXYO2)  &
    + rct(60,k) * xnew(C5DICARBO2)  &
    + rct(60,k) * xnew(ISRO2)  &
    + rct(47,k) * xnew(MACRO2)
  
  !-> HO2 cont.
    L = L +   rct(76,k) * xnew(TERPO2)  &
    + rct(100,k)
  
  xnew(HO2) = (xold(HO2) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3O2

  P = rct(27,k) * xnew(OH)  &
    + 0.6* rct(34,k) * xnew(CH3OOH) * xnew(OH)  &
    + 0.05* rct(40,k) * xnew(CH3CHO) * xnew(OH)  &
    + rct(42,k) * xnew(CH3CO3) * xnew(NO)  &
    + 4.0e-12 * xnew(CH3CO3) * xnew(NO3)  &
    + 0.454* rct(45,k) * xnew(CH3CO3) * xnew(HO2)  &
    + 0.7* 1.0e-11 * xnew(RO2POOL) * xnew(CH3CO3)  &
    + 0.3* rct(53,k) * xnew(C3H6) * xnew(O3)  &
    + 0.08* rct(65,k) * xnew(C5H8) * xnew(O3)  &
    + rcphot(IDCH3CHO,k) * xnew(CH3CHO)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3CO3H)
  
  L = 1.2e-12 * xnew(NO3)  &
    + rct(29,k) * xnew(NO)  &
    + rct(30,k) * xnew(RO2POOL)  &
    + rct(31,k) * xnew(RO2POOL)  &
    + rct(32,k) * xnew(HO2)
  
  xnew(CH3O2) = (xold(CH3O2) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5O2

  P = rct(37,k) * xnew(C2H6) * xnew(OH)  &
    + 0.31* 1.17e-11 * xnew(C2H5OOH) * xnew(OH)  &
    + 0.264* rct(48,k) * xnew(C4H9O2) * xnew(NO)  &
    + 0.284* 2.3e-12 * xnew(C4H9O2) * xnew(NO3)  &
    + rct(55,k) * xnew(C3H6) * xnew(NO3)  &
    + rcphot(IDMEK,k) * xnew(MEK)  &
    + 0.284* rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + 0.348* rcphot(IDMEK,k) *0.78 * xnew(SC4H9NO3)
  
  L = rct(38,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + rct(39,k) * xnew(HO2)
  
  xnew(C2H5O2) = (xold(C2H5O2) + dt2 * P) / (1.0 + dt2 * L)


!-> C4H9O2

  P = rct(46,k) * xnew(NC4H10) * xnew(OH)  &
    + 0.153* 2.4e-11 * xnew(BURO2H) * xnew(OH)
  
  L = rct(47,k) * xnew(HO2)  &
    + rct(48,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)
  
  xnew(C4H9O2) = (xold(C4H9O2) + dt2 * P) / (1.0 + dt2 * L)


!-> ISRO2

  P = rct(64,k) * xnew(C5H8) * xnew(OH)
  
  L = rct(48,k) * xnew(NO)  &
    + rct(60,k) * xnew(HO2)  &
    + 4.0e-12 * xnew(RO2POOL)  &
    + rct(73,k)  &
    + rct(74,k)
  
  xnew(ISRO2) = (xold(ISRO2) + dt2 * P) / (1.0 + dt2 * L)


!-> ETRO2

  P = 0.05* 3.22e-12 * xnew(C2H5OH) * xnew(OH)  &
    + rct(49,k) * xnew(C2H4) * xnew(OH)  &
    + 0.21* 1.7475e-11 * xnew(ETRO2H) * xnew(OH)
  
  L = rct(48,k) * xnew(NO)  &
    + rct(50,k) * xnew(HO2)
  
  xnew(ETRO2) = (xold(ETRO2) + dt2 * P) / (1.0 + dt2 * L)


!-> PRRO2

  P = rct(52,k) * xnew(C3H6) * xnew(OH)
  
  L = rct(48,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + rct(54,k) * xnew(HO2)
  
  xnew(PRRO2) = (xold(PRRO2) + dt2 * P) / (1.0 + dt2 * L)


!-> OXYO2

  P = 0.259* rct(56,k) * xnew(BENZENE) * xnew(OH)  &
    + 0.564* rct(57,k) * xnew(TOLUENE) * xnew(OH)  &
    + 0.55* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + 7.62e-11 * xnew(OXYO2H) * xnew(OH)
  
  L = rct(48,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + rct(58,k) * xnew(HO2)
  
  xnew(OXYO2) = (xold(OXYO2) + dt2 * P) / (1.0 + dt2 * L)


!-> MEKO2

  P = 1.1e-12 * xnew(MEK) * xnew(OH)  &
    + 0.174* 2.1e-11 * xnew(MEKO2H) * xnew(OH)
  
  L = rct(48,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + rct(47,k) * xnew(HO2)  &
    + 1.04e-12 * xnew(RO2POOL)
  
  xnew(MEKO2) = (xold(MEKO2) + dt2 * P) / (1.0 + dt2 * L)


!-> C5DICARBO2

  P = 6.2e-11 * xnew(C5DICARB) * xnew(OH)  &
    + rct(61,k) * xnew(C5DICAROOH) * xnew(OH)  &
    + rcphot(IDNO2,k) *0.2 * xnew(C5DICARB)
  
  L = rct(48,k) * xnew(NO)  &
    + rct(59,k) * xnew(NO3)  &
    + rct(60,k) * xnew(HO2)
  
  xnew(C5DICARBO2) = (xold(C5DICARBO2) + dt2 * P) / (1.0 + dt2 * L)


!-> MACRO2

  P = 0.1* rct(65,k) * xnew(C5H8) * xnew(O3)  &
    + 9.13e-12 * xnew(IEPOX) * xnew(OH)  &
    + rct(67,k) * xnew(MACR) * xnew(OH)  &
    + rct(68,k) * xnew(MACR) * xnew(OH)  &
    + rct(44,k) * xnew(MPAN)
  
  L = rct(48,k) * xnew(NO)  &
    + rct(47,k) * xnew(HO2)  &
    + rct(72,k) * xnew(NO2)
  
  xnew(MACRO2) = (xold(MACRO2) + dt2 * P) / (1.0 + dt2 * L)


!-> TERPO2

  P = rct(75,k) * xnew(APINENE) * xnew(OH)  &
    + rct(77,k) * xnew(TERPOOH) * xnew(OH)  &
    + rct(79,k) * xnew(APINENE) * xnew(NO3)
  
  L = rct(48,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + rct(76,k) * xnew(HO2)
  
  xnew(TERPO2) = (xold(TERPO2) + dt2 * P) / (1.0 + dt2 * L)


!-> RO2POOL

  ! P = 0.0
  
  ! L = 0.0
  
  xnew(RO2POOL) = sum(xnew(RO2_GROUP))


!-> CH3CO3

  P = 0.95* rct(40,k) * xnew(CH3CHO) * xnew(OH)  &
    + rct(41,k) * xnew(CH3CHO) * xnew(NO3)  &
    + rct(44,k) * xnew(PAN)  &
    + 1.1e-11 * xnew(CH3CO3H) * xnew(OH)  &
    + rct(48,k) * xnew(MEKO2) * xnew(NO)  &
    + 2.3e-12 * xnew(MEKO2) * xnew(NO3)  &
    + 0.75* 1.04e-12 * xnew(RO2POOL) * xnew(MEKO2)  &
    + 2.3e-13 * xnew(BIACET) * xnew(OH)  &
    + rct(63,k) * xnew(MGLYOX) * xnew(OH)  &
    + 0.1* rct(65,k) * xnew(C5H8) * xnew(O3)  &
    + 0.1* rct(69,k) * xnew(MACR) * xnew(O3)  &
    + 0.1* rct(70,k) * xnew(MACR) * xnew(O3)  &
    + 0.25* rct(48,k) * xnew(MACRO2) * xnew(NO)  &
    + rcphot(IDMEK,k) * xnew(MEK)  &
    + rcphot(IDMEK,k) * xnew(ACETOL)  &
    + rcphot(IDRCOCHO,k) * xnew(MGLYOX)  &
    + 2* rcphot(IDCH3COY,k) * xnew(BIACET)  &
    + rcphot(IDHCHO_H2,k) *0.0474 * xnew(MACR)  &
    + rcphot(IDCH3O2H,k) * xnew(MEKO2H)  &
    + rcphot(IDMEK,k) * xnew(MEKO2H)
  
  !-> CH3CO3 cont.
    P = P +   rcphot(IDPAN,k) * xnew(PAN)
  
  L = rct(42,k) * xnew(NO)  &
    + 4.0e-12 * xnew(NO3)  &
    + rct(43,k) * xnew(NO2)  &
    + rct(45,k) * xnew(HO2)  &
    + 1.0e-11 * xnew(RO2POOL)
  
  xnew(CH3CO3) = (xold(CH3CO3) + dt2 * P) / (1.0 + dt2 * L)


!-> O3

  P = rct(1,k) * xnew(OP)  &
    + 0.162* rct(45,k) * xnew(CH3CO3) * xnew(HO2)
  
  L = rct(11,k) * xnew(NO)  &
    + rct(12,k) * xnew(NO2)  &
    + rct(13,k) * xnew(OH)  &
    + rct(14,k) * xnew(HO2)  &
    + rct(51,k) * xnew(C2H4)  &
    + rct(53,k) * xnew(C3H6)  &
    + rct(65,k) * xnew(C5H8)  &
    + 9.2e-17 * xnew(ISON)  &
    + rct(69,k) * xnew(MACR)  &
    + rct(70,k) * xnew(MACR)  &
    + rct(78,k) * xnew(APINENE)  &
    + rcphot(IDO3_O1D,k)  &
    + rcphot(IDO3_O3P,k)  &
    + AQRCK(ICLRC2,K) * xnew(SO2)
  
  xnew(O3) = (xold(O3) + dt2 * P) / (1.0 + dt2 * L)


!-> NO

  P = rct(17,k) * xnew(NO2) * xnew(NO3)  &
    + rcphot(IDNO2,k) * xnew(NO2)  &
    + 0.127* rcphot(IDNO3,k) * xnew(NO3)  &
    + rcphot(IDHONO,k) * xnew(HONO)  &
    + rcemis(NO,k)
  
  L = rct(10,k) * xnew(OP)  &
    + rct(11,k) * xnew(O3)  &
    + rct(15,k) * xnew(NO3)  &
    + rct(16,k) * xnew(HO2)  &
    + rct(22,k) * xnew(OH)  &
    + rct(29,k) * xnew(CH3O2)  &
    + rct(38,k) * xnew(C2H5O2)  &
    + rct(42,k) * xnew(CH3CO3)  &
    + rct(48,k) * xnew(C4H9O2)  &
    + rct(48,k) * xnew(MEKO2)  &
    + rct(48,k) * xnew(ETRO2)  &
    + rct(48,k) * xnew(PRRO2)  &
    + rct(48,k) * xnew(OXYO2)  &
    + rct(48,k) * xnew(C5DICARBO2)  &
    + rct(48,k) * xnew(ISRO2)  &
    + rct(48,k) * xnew(MACRO2)  &
    + rct(48,k) * xnew(TERPO2)
  
  xnew(NO) = (xold(NO) + dt2 * P) / (1.0 + dt2 * L)


!-> NO2

  P = rct(10,k) * xnew(OP) * xnew(NO)  &
    + rct(11,k) * xnew(O3) * xnew(NO)  &
    + 2* rct(15,k) * xnew(NO) * xnew(NO3)  &
    + rct(16,k) * xnew(NO) * xnew(HO2)  &
    + rct(17,k) * xnew(NO2) * xnew(NO3)  &
    + rct(19,k) * xnew(N2O5)  &
    + rct(23,k) * xnew(OH) * xnew(HONO)  &
    + 2.0e-11 * xnew(OH) * xnew(NO3)  &
    + 4.0e-12 * xnew(HO2) * xnew(NO3)  &
    + rct(25,k) * xnew(HO2NO2)  &
    + rct(26,k) * xnew(OH) * xnew(HO2NO2)  &
    + 1.2e-12 * xnew(CH3O2) * xnew(NO3)  &
    + rct(29,k) * xnew(CH3O2) * xnew(NO)  &
    + rct(38,k) * xnew(C2H5O2) * xnew(NO)  &
    + 2.3e-12 * xnew(C2H5O2) * xnew(NO3)  &
    + rct(42,k) * xnew(CH3CO3) * xnew(NO)  &
    + 4.0e-12 * xnew(CH3CO3) * xnew(NO3)  &
    + rct(44,k) * xnew(PAN)  &
    + 0.917* rct(48,k) * xnew(C4H9O2) * xnew(NO)  &
    + 2.3e-12 * xnew(C4H9O2) * xnew(NO3)
  
  !-> NO2 cont.
    P = P +   8.6e-13 * xnew(SC4H9NO3) * xnew(OH)  &
    + rct(48,k) * xnew(MEKO2) * xnew(NO)  &
    + 2.3e-12 * xnew(MEKO2) * xnew(NO3)  &
    + rct(48,k) * xnew(ETRO2) * xnew(NO)  &
    + rct(48,k) * xnew(PRRO2) * xnew(NO)  &
    + 2.3e-12 * xnew(PRRO2) * xnew(NO3)  &
    + rct(55,k) * xnew(C3H6) * xnew(NO3)  &
    + rct(48,k) * xnew(OXYO2) * xnew(NO)  &
    + 2.3e-12 * xnew(OXYO2) * xnew(NO3)  &
    + rct(48,k) * xnew(C5DICARBO2) * xnew(NO)  &
    + rct(59,k) * xnew(C5DICARBO2) * xnew(NO3)  &
    + 0.9* rct(48,k) * xnew(ISRO2) * xnew(NO)  &
    + 0.22* 3.34e-11 * xnew(ISON) * xnew(OH)  &
    + 0.33* 9.2e-17 * xnew(ISON) * xnew(O3)  &
    + rct(48,k) * xnew(MACRO2) * xnew(NO)  &
    + 2.9e-11 * xnew(MPAN) * xnew(OH)  &
    + 1.55e-12 * xnew(NALD) * xnew(OH)  &
    + rct(44,k) * xnew(MPAN)  &
    + rct(48,k) * xnew(TERPO2) * xnew(NO)  &
    + 2.3e-12 * xnew(TERPO2) * xnew(NO3)
  
  !-> NO2 cont.
    P = P +   rct(79,k) * xnew(APINENE) * xnew(NO3)  &
    + 0.873* rcphot(IDNO3,k) * xnew(NO3)  &
    + rcphot(IDHNO3,k) * xnew(HNO3)  &
    + 0.667* rcphot(IDHO2NO2,k) * xnew(HO2NO2)  &
    + rcphot(IDCH3CHO,k) * xnew(NALD)  &
    + rcphot(IDMEK,k) *0.78 * xnew(SC4H9NO3)  &
    + rcphot(IDMEK,k) *0.78 * xnew(ISON)  &
    + 0.7* rcphot(IDPAN,k) * xnew(PAN)  &
    + rcemis(NO2,k)
  
  L = rct(12,k) * xnew(O3)  &
    + rct(17,k) * xnew(NO3)  &
    + rct(18,k) * xnew(NO3)  &
    + rct(20,k) * xnew(OH)  &
    + rct(24,k) * xnew(HO2)  &
    + rct(43,k) * xnew(CH3CO3)  &
    + rct(72,k) * xnew(MACRO2)  &
    + rcphot(IDNO2,k)
  
  xnew(NO2) = (xold(NO2) + dt2 * P) / (1.0 + dt2 * L)


!-> NO3

  P = rct(12,k) * xnew(O3) * xnew(NO2)  &
    + rct(19,k) * xnew(N2O5)  &
    + rct(21,k) * xnew(OH) * xnew(HNO3)  &
    + 0.333* rcphot(IDHO2NO2,k) * xnew(HO2NO2)  &
    + 0.3* rcphot(IDPAN,k) * xnew(PAN)
  
  L = rct(15,k) * xnew(NO)  &
    + rct(17,k) * xnew(NO2)  &
    + rct(18,k) * xnew(NO2)  &
    + 2.0e-11 * xnew(OH)  &
    + 4.0e-12 * xnew(HO2)  &
    + 1.2e-12 * xnew(CH3O2)  &
    + rct(36,k) * xnew(HCHO)  &
    + 2.3e-12 * xnew(C2H5O2)  &
    + rct(41,k) * xnew(CH3CHO)  &
    + 4.0e-12 * xnew(CH3CO3)  &
    + 2.3e-12 * xnew(C4H9O2)  &
    + 2.3e-12 * xnew(MEKO2)  &
    + 2.3e-12 * xnew(PRRO2)  &
    + rct(55,k) * xnew(C3H6)  &
    + 2.3e-12 * xnew(OXYO2)  &
    + rct(59,k) * xnew(C5DICARBO2)  &
    + rct(66,k) * xnew(C5H8)  &
    + 2.3e-12 * xnew(TERPO2)  &
    + rct(79,k) * xnew(APINENE)  &
    + rcphot(IDNO3,k)
  
  !-> NO3 cont.
    L = L +   rct(97,k)
  
  xnew(NO3) = (xold(NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> N2O5

  P = rct(18,k) * xnew(NO2) * xnew(NO3)
  
  L = rct(19,k)  &
    + rct(96,k)
  
  xnew(N2O5) = (xold(N2O5) + dt2 * P) / (1.0 + dt2 * L)


!-> H2

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for H2! xnew(H2) = max(0.0, xold(H2))


!-> H2O2

  P = rct(8,k) * xnew(HO2) * xnew(HO2)  &
    + rct(9,k) * xnew(HO2) * xnew(HO2)  &
    + 0.1575* rct(51,k) * xnew(C2H4) * xnew(O3)  &
    + 0.169* rct(53,k) * xnew(C3H6) * xnew(O3)  &
    + 0.09* rct(65,k) * xnew(C5H8) * xnew(O3)  &
    + 0.5* rct(100,k) * xnew(HO2)
  
  L = rct(6,k) * xnew(OH)  &
    + rcphot(IDH2O2,k)  &
    + AQRCK(ICLRC1,K) * xnew(SO2)  &
    + 1.0e-5
  
  xnew(H2O2) = (xold(H2O2) + dt2 * P) / (1.0 + dt2 * L)


!-> HONO

  P = rct(22,k) * xnew(OH) * xnew(NO)
  
  L = rct(23,k) * xnew(OH)  &
    + rcphot(IDHONO,k)
  
  xnew(HONO) = (xold(HONO) + dt2 * P) / (1.0 + dt2 * L)


!-> HNO3

  P = rct(20,k) * xnew(NO2) * xnew(OH)  &
    + rct(36,k) * xnew(HCHO) * xnew(NO3)  &
    + rct(41,k) * xnew(CH3CHO) * xnew(NO3)  &
    + 2* rct(96,k) * xnew(N2O5)  &
    + rct(97,k) * xnew(NO3)  &
    + rct(20,k) * xnew(shipNOx) * xnew(OH)  &
    + 3.2e-5 * xnew(shipNOx)
  
  L = rct(21,k) * xnew(OH)  &
    + rcphot(IDHNO3,k)  &
    + rct(98,k)  &
    + rct(99,k)
  
  xnew(HNO3) = (xold(HNO3) + dt2 * P) / (1.0 + dt2 * L)


!-> HO2NO2

  P = rct(24,k) * xnew(HO2) * xnew(NO2)
  
  L = rct(25,k)  &
    + rct(26,k) * xnew(OH)  &
    + rcphot(IDHO2NO2,k)
  
  xnew(HO2NO2) = (xold(HO2NO2) + dt2 * P) / (1.0 + dt2 * L)


!-> CO

  P = rct(35,k) * xnew(HCHO) * xnew(OH)  &
    + rct(36,k) * xnew(HCHO) * xnew(NO3)  &
    + 0.05* rct(40,k) * xnew(CH3CHO) * xnew(OH)  &
    + 2.3e-13 * xnew(BIACET) * xnew(OH)  &
    + 0.35* rct(51,k) * xnew(C2H4) * xnew(O3)  &
    + 0.483* rct(53,k) * xnew(C3H6) * xnew(O3)  &
    + 1.416* rct(62,k) * xnew(GLYOX) * xnew(OH)  &
    + rct(63,k) * xnew(MGLYOX) * xnew(OH)  &
    + 0.14* rct(65,k) * xnew(C5H8) * xnew(O3)  &
    + 0.22* rct(69,k) * xnew(MACR) * xnew(O3)  &
    + 0.22* rct(70,k) * xnew(MACR) * xnew(O3)  &
    + 0.25* rct(48,k) * xnew(MACRO2) * xnew(NO)  &
    + 3.0e-11 * xnew(MACROOH) * xnew(OH)  &
    + 1.55e-12 * xnew(NALD) * xnew(OH)  &
    + 7.61e-11 * xnew(HPALD) * xnew(OH)  &
    + rcphot(IDHCHO_H,k) * xnew(HCHO)  &
    + rcphot(IDHCHO_H2,k) * xnew(HCHO)  &
    + rcphot(IDCH3CHO,k) * xnew(CH3CHO)  &
    + rcphot(IDCH3CHO,k) * xnew(NALD)  &
    + rcphot(IDRCOCHO,k) * xnew(MGLYOX)
  
  !-> CO cont.
    P = P +   rcphot(IDHCHO_H2,k) *0.0474 * xnew(MACR)  &
    + rcphot(IDHCHO_H2,k) *24. * xnew(PACALD)  &
    + 0.5* rcphot(IDCH3O2H,k) * xnew(MACROOH)  &
    + 2* rcphot(IDGLYOXA,k) * xnew(GLYOX)  &
    + 2* rcphot(IDGLYOXB,k) * xnew(GLYOX)  &
    + rcphot(IDGLYOXC,k) * xnew(GLYOX)  &
    + rcemis(CO,k)
  
  L = rct(28,k) * xnew(OH)
  
  xnew(CO) = (xold(CO) + dt2 * P) / (1.0 + dt2 * L)


!-> CH4

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for CH4! xnew(CH4) = max(0.0, xold(CH4))


!-> C2H6

  P = 0.073* 1.04e-12 * xnew(RO2POOL) * xnew(MEKO2)  &
    + 0.0603* 2.1e-11 * xnew(MEKO2H) * xnew(OH)  &
    + 0.19* rct(48,k) * xnew(TERPO2) * xnew(NO)  &
    + 0.19* 2.3e-12 * xnew(TERPO2) * xnew(NO3)  &
    + rcemis(C2H6,k)
  
  L = rct(37,k) * xnew(OH)
  
  xnew(C2H6) = (xold(C2H6) + dt2 * P) / (1.0 + dt2 * L)


!-> NC4H10

  P = rcemis(NC4H10,k)
  
  L = rct(46,k) * xnew(OH)
  
  xnew(NC4H10) = (xold(NC4H10) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H4

  P = rcemis(C2H4,k)
  
  L = rct(49,k) * xnew(OH)  &
    + rct(51,k) * xnew(O3)
  
  xnew(C2H4) = (xold(C2H4) + dt2 * P) / (1.0 + dt2 * L)


!-> C3H6

  P = rcemis(C3H6,k)
  
  L = rct(52,k) * xnew(OH)  &
    + rct(53,k) * xnew(O3)  &
    + rct(55,k) * xnew(NO3)
  
  xnew(C3H6) = (xold(C3H6) + dt2 * P) / (1.0 + dt2 * L)


!-> BENZENE

  P = rcemis(BENZENE,k)
  
  L = rct(56,k) * xnew(OH)
  
  xnew(BENZENE) = (xold(BENZENE) + dt2 * P) / (1.0 + dt2 * L)


!-> TOLUENE

  P = rcemis(TOLUENE,k)
  
  L = rct(57,k) * xnew(OH)
  
  xnew(TOLUENE) = (xold(TOLUENE) + dt2 * P) / (1.0 + dt2 * L)


!-> OXYL

  P = rcemis(OXYL,k)
  
  L = 1.36e-11 * xnew(OH)
  
  xnew(OXYL) = (xold(OXYL) + dt2 * P) / (1.0 + dt2 * L)


!-> C5H8

  P = rcemis(C5H8,k)  &
    + RCBIO(1,K)
  
  L = rct(64,k) * xnew(OH)  &
    + rct(65,k) * xnew(O3)  &
    + rct(66,k) * xnew(NO3)
  
  xnew(C5H8) = (xold(C5H8) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3OH

  P = 0.5* rct(31,k) * xnew(RO2POOL) * xnew(CH3O2)  &
    + 0.05* rct(53,k) * xnew(C3H6) * xnew(O3)  &
    + rcemis(CH3OH,k)
  
  L = rct(33,k) * xnew(OH)
  
  xnew(CH3OH) = (xold(CH3OH) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5OH

  P = rcemis(C2H5OH,k)
  
  L = 3.22e-12 * xnew(OH)
  
  xnew(C2H5OH) = (xold(C2H5OH) + dt2 * P) / (1.0 + dt2 * L)


!-> HCHO

  P = 1.2e-12 * xnew(CH3O2) * xnew(NO3)  &
    + rct(29,k) * xnew(CH3O2) * xnew(NO)  &
    + rct(30,k) * xnew(RO2POOL) * xnew(CH3O2)  &
    + 0.5* rct(31,k) * xnew(RO2POOL) * xnew(CH3O2)  &
    + 0.1* rct(32,k) * xnew(CH3O2) * xnew(HO2)  &
    + rct(33,k) * xnew(CH3OH) * xnew(OH)  &
    + 0.4* rct(34,k) * xnew(CH3OOH) * xnew(OH)  &
    + 0.1095* 1.04e-12 * xnew(RO2POOL) * xnew(MEKO2)  &
    + 2.3e-13 * xnew(BIACET) * xnew(OH)  &
    + 2* rct(48,k) * xnew(ETRO2) * xnew(NO)  &
    + 1.1575* rct(51,k) * xnew(C2H4) * xnew(O3)  &
    + 0.68* rct(53,k) * xnew(C3H6) * xnew(O3)  &
    + rct(48,k) * xnew(PRRO2) * xnew(NO)  &
    + 2.3e-12 * xnew(PRRO2) * xnew(NO3)  &
    + rct(55,k) * xnew(C3H6) * xnew(NO3)  &
    + 0.58* rct(65,k) * xnew(C5H8) * xnew(O3)  &
    + 0.9* rct(48,k) * xnew(ISRO2) * xnew(NO)  &
    + 4.0e-12 * xnew(RO2POOL) * xnew(ISRO2)  &
    + 0.22* 3.34e-11 * xnew(ISON) * xnew(OH)  &
    + 0.33* 9.2e-17 * xnew(ISON) * xnew(O3)
  
  !-> HCHO cont.
    P = P +   0.75* rct(48,k) * xnew(MACRO2) * xnew(NO)  &
    + 1.55e-12 * xnew(NALD) * xnew(OH)  &
    + rct(74,k) * xnew(ISRO2)  &
    + 7.61e-11 * xnew(HPALD) * xnew(OH)  &
    + rcphot(IDCH3CHO,k) * xnew(NALD)  &
    + rcphot(IDMEK,k) * xnew(ACETOL)  &
    + rcphot(IDHCHO_H2,k) *0.0474 * xnew(MACR)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3OOH)  &
    + 1.56* rcphot(IDCH3O2H,k) * xnew(ETRO2H)  &
    + 0.5* rcphot(IDCH3O2H,k) * xnew(MACROOH)  &
    + rcphot(IDCH3O2H,k) * xnew(ISRO2H)  &
    + rcphot(IDMEK,k) *0.78 * xnew(ISON)  &
    + rcphot(IDGLYOXC,k) * xnew(GLYOX)  &
    + rcemis(HCHO,k)
  
  L = rct(35,k) * xnew(OH)  &
    + rct(36,k) * xnew(NO3)  &
    + rcphot(IDHCHO_H,k)  &
    + rcphot(IDHCHO_H2,k)
  
  xnew(HCHO) = (xold(HCHO) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3CHO

  P = rct(38,k) * xnew(C2H5O2) * xnew(NO)  &
    + 2.3e-12 * xnew(C2H5O2) * xnew(NO3)  &
    + 0.69* 1.17e-11 * xnew(C2H5OOH) * xnew(OH)  &
    + 0.95* 3.22e-12 * xnew(C2H5OH) * xnew(OH)  &
    + 0.264* rct(48,k) * xnew(C4H9O2) * xnew(NO)  &
    + 0.284* 2.3e-12 * xnew(C4H9O2) * xnew(NO3)  &
    + rct(48,k) * xnew(MEKO2) * xnew(NO)  &
    + 2.3e-12 * xnew(MEKO2) * xnew(NO3)  &
    + 0.6405* 1.04e-12 * xnew(RO2POOL) * xnew(MEKO2)  &
    + 0.79* 1.7475e-11 * xnew(ETRO2H) * xnew(OH)  &
    + 0.489* rct(53,k) * xnew(C3H6) * xnew(O3)  &
    + rct(48,k) * xnew(PRRO2) * xnew(NO)  &
    + 2.3e-12 * xnew(PRRO2) * xnew(NO3)  &
    + rcphot(IDCH3O2H,k) *0.5 * xnew(C2H5OOH)  &
    + 0.284* rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + 0.22* rcphot(IDCH3O2H,k) * xnew(ETRO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(MEKO2H)  &
    + rcphot(IDMEK,k) * xnew(MEKO2H)  &
    + 0.348* rcphot(IDMEK,k) *0.78 * xnew(SC4H9NO3)  &
    + rcemis(CH3CHO,k)
  
  L = rct(40,k) * xnew(OH)  &
    + rct(41,k) * xnew(NO3)  &
    + rcphot(IDCH3CHO,k)
  
  xnew(CH3CHO) = (xold(CH3CHO) + dt2 * P) / (1.0 + dt2 * L)


!-> MACR

  P = 0.65* rct(65,k) * xnew(C5H8) * xnew(O3)  &
    + 0.9* rct(48,k) * xnew(ISRO2) * xnew(NO)  &
    + 4.0e-12 * xnew(RO2POOL) * xnew(ISRO2)  &
    + 0.1* 8.94e-11 * xnew(ISRO2H) * xnew(OH)  &
    + 0.22* 3.34e-11 * xnew(ISON) * xnew(OH)  &
    + 0.33* 9.2e-17 * xnew(ISON) * xnew(O3)  &
    + rct(74,k) * xnew(ISRO2)  &
    + 2* rct(48,k) * xnew(TERPO2) * xnew(NO)  &
    + 2* 2.3e-12 * xnew(TERPO2) * xnew(NO3)  &
    + 2* rct(78,k) * xnew(APINENE) * xnew(O3)  &
    + rcphot(IDCH3O2H,k) * xnew(ISRO2H)  &
    + rcphot(IDMEK,k) *0.78 * xnew(ISON)
  
  L = rct(67,k) * xnew(OH)  &
    + rct(68,k) * xnew(OH)  &
    + rct(69,k) * xnew(O3)  &
    + rct(70,k) * xnew(O3)  &
    + rcphot(IDHCHO_H2,k) *0.0474
  
  xnew(MACR) = (xold(MACR) + dt2 * P) / (1.0 + dt2 * L)


!-> MEK

  P = 0.847* 2.4e-11 * xnew(BURO2H) * xnew(OH)  &
    + 0.654* rct(48,k) * xnew(C4H9O2) * xnew(NO)  &
    + 0.716* 2.3e-12 * xnew(C4H9O2) * xnew(NO3)  &
    + 8.6e-13 * xnew(SC4H9NO3) * xnew(OH)  &
    + 0.716* rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + 0.652* rcphot(IDMEK,k) *0.78 * xnew(SC4H9NO3)  &
    + rcemis(MEK,k)
  
  L = 1.1e-12 * xnew(OH)  &
    + rcphot(IDMEK,k)
  
  xnew(MEK) = (xold(MEK) + dt2 * P) / (1.0 + dt2 * L)


!-> ACETOL

  P = rct(54,k) * xnew(PRRO2) * xnew(HO2)  &
    + 0.78* 3.34e-11 * xnew(ISON) * xnew(OH)  &
    + 0.25* rct(48,k) * xnew(MACRO2) * xnew(NO)  &
    + 2.9e-11 * xnew(MPAN) * xnew(OH)  &
    + 3.0e-11 * xnew(MACROOH) * xnew(OH)  &
    + 0.5* rcphot(IDCH3O2H,k) * xnew(MACROOH)
  
  L = rct(71,k) * xnew(OH)  &
    + rcphot(IDMEK,k)
  
  xnew(ACETOL) = (xold(ACETOL) + dt2 * P) / (1.0 + dt2 * L)


!-> GLYOX

  P = rct(48,k) * xnew(C5DICARBO2) * xnew(NO)  &
    + rct(59,k) * xnew(C5DICARBO2) * xnew(NO3)  &
    + rcphot(IDHCHO_H,k) *0.9 * xnew(C5DICAROOH)  &
    + rcemis(GLYOX,k)
  
  L = rct(62,k) * xnew(OH)  &
    + rcphot(IDGLYOXA,k)  &
    + rcphot(IDGLYOXB,k)  &
    + rcphot(IDGLYOXC,k)
  
  xnew(GLYOX) = (xold(GLYOX) + dt2 * P) / (1.0 + dt2 * L)


!-> MGLYOX

  P = 0.0365* 1.04e-12 * xnew(RO2POOL) * xnew(MEKO2)  &
    + 0.1206* 2.1e-11 * xnew(MEKO2H) * xnew(OH)  &
    + rct(48,k) * xnew(OXYO2) * xnew(NO)  &
    + 2.3e-12 * xnew(OXYO2) * xnew(NO3)  &
    + rct(48,k) * xnew(C5DICARBO2) * xnew(NO)  &
    + rct(59,k) * xnew(C5DICARBO2) * xnew(NO3)  &
    + 0.9* rct(69,k) * xnew(MACR) * xnew(O3)  &
    + 0.9* rct(70,k) * xnew(MACR) * xnew(O3)  &
    + 0.5* rct(48,k) * xnew(MACRO2) * xnew(NO)  &
    + rct(71,k) * xnew(ACETOL) * xnew(OH)  &
    + 7.61e-11 * xnew(HPALD) * xnew(OH)  &
    + rcphot(IDHCHO_H2,k) *24. * xnew(PACALD)  &
    + rcphot(IDHCHO_H,k) *0.9 * xnew(C5DICAROOH)  &
    + 0.5* rcphot(IDCH3O2H,k) * xnew(MACROOH)  &
    + rcphot(IDCH3O2H,k) * xnew(OXYO2H)  &
    + rcemis(MGLYOX,k)
  
  L = rct(63,k) * xnew(OH)  &
    + rcphot(IDRCOCHO,k)
  
  xnew(MGLYOX) = (xold(MGLYOX) + dt2 * P) / (1.0 + dt2 * L)


!-> BIACET

  P = 0.2135* 1.04e-12 * xnew(RO2POOL) * xnew(MEKO2)  &
    + 0.7054* 2.1e-11 * xnew(MEKO2H) * xnew(OH)
  
  L = 2.3e-13 * xnew(OH)  &
    + rcphot(IDCH3COY,k)
  
  xnew(BIACET) = (xold(BIACET) + dt2 * P) / (1.0 + dt2 * L)


!-> C5DICARB

  P = rct(48,k) * xnew(OXYO2) * xnew(NO)  &
    + 2.3e-12 * xnew(OXYO2) * xnew(NO3)  &
    + rcphot(IDCH3O2H,k) * xnew(OXYO2H)
  
  L = 6.2e-11 * xnew(OH)  &
    + rcphot(IDNO2,k) *0.2
  
  xnew(C5DICARB) = (xold(C5DICARB) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3OOH

  P = 0.9* rct(32,k) * xnew(CH3O2) * xnew(HO2)
  
  L = rct(34,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)  &
    + 1.0e-5
  
  xnew(CH3OOH) = (xold(CH3OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5OOH

  P = rct(39,k) * xnew(C2H5O2) * xnew(HO2)
  
  L = 1.17e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k) *0.5
  
  xnew(C2H5OOH) = (xold(C2H5OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> BURO2H

  P = rct(47,k) * xnew(C4H9O2) * xnew(HO2)
  
  L = 2.4e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(BURO2H) = (xold(BURO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> ETRO2H

  P = rct(50,k) * xnew(ETRO2) * xnew(HO2)
  
  L = 1.7475e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(ETRO2H) = (xold(ETRO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> MEKO2H

  P = rct(47,k) * xnew(MEKO2) * xnew(HO2)
  
  L = 2.1e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)  &
    + rcphot(IDMEK,k)
  
  xnew(MEKO2H) = (xold(MEKO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> ISRO2H

  P = rct(60,k) * xnew(ISRO2) * xnew(HO2)
  
  L = 8.94e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(ISRO2H) = (xold(ISRO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> C5DICAROOH

  P = rct(60,k) * xnew(C5DICARBO2) * xnew(HO2)
  
  L = rct(61,k) * xnew(OH)  &
    + rcphot(IDHCHO_H,k) *0.9
  
  xnew(C5DICAROOH) = (xold(C5DICAROOH) + dt2 * P) / (1.0 + dt2 * L)


!-> HPALD

  P = rct(73,k) * xnew(ISRO2)
  
  L = 7.61e-11 * xnew(OH)  &
    + rcphot(IDHCHO_H2,k) *12.
  
  xnew(HPALD) = (xold(HPALD) + dt2 * P) / (1.0 + dt2 * L)


!-> MACROOH

  P = rct(47,k) * xnew(MACRO2) * xnew(HO2)
  
  L = 3.0e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(MACROOH) = (xold(MACROOH) + dt2 * P) / (1.0 + dt2 * L)


!-> OXYO2H

  P = rct(58,k) * xnew(OXYO2) * xnew(HO2)
  
  L = 7.62e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(OXYO2H) = (xold(OXYO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3CO3H

  P = 0.384* rct(45,k) * xnew(CH3CO3) * xnew(HO2)
  
  L = 1.1e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(CH3CO3H) = (xold(CH3CO3H) + dt2 * P) / (1.0 + dt2 * L)


!-> PACALD

  P = rcphot(IDHCHO_H2,k) *12. * xnew(HPALD)
  
  L = rcphot(IDHCHO_H2,k) *24.
  
  xnew(PACALD) = (xold(PACALD) + dt2 * P) / (1.0 + dt2 * L)


!-> IEPOX

  P = 0.9* 8.94e-11 * xnew(ISRO2H) * xnew(OH)
  
  L = 9.13e-12 * xnew(OH)
  
  xnew(IEPOX) = (xold(IEPOX) + dt2 * P) / (1.0 + dt2 * L)


!-> SC4H9NO3

  P = 0.083* rct(48,k) * xnew(C4H9O2) * xnew(NO)
  
  L = 8.6e-13 * xnew(OH)  &
    + rcphot(IDMEK,k) *0.78
  
  xnew(SC4H9NO3) = (xold(SC4H9NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> NALD

  P = 0.78* 3.34e-11 * xnew(ISON) * xnew(OH)  &
    + 0.67* 9.2e-17 * xnew(ISON) * xnew(O3)
  
  L = 1.55e-12 * xnew(OH)  &
    + rcphot(IDCH3CHO,k)
  
  xnew(NALD) = (xold(NALD) + dt2 * P) / (1.0 + dt2 * L)


!-> ISON

  P = rct(66,k) * xnew(C5H8) * xnew(NO3)  &
    + 0.1* rct(48,k) * xnew(ISRO2) * xnew(NO)
  
  L = 3.34e-11 * xnew(OH)  &
    + 9.2e-17 * xnew(O3)  &
    + rcphot(IDMEK,k) *0.78
  
  xnew(ISON) = (xold(ISON) + dt2 * P) / (1.0 + dt2 * L)


!-> PAN

  P = rct(43,k) * xnew(CH3CO3) * xnew(NO2)
  
  L = rct(44,k)  &
    + rcphot(IDPAN,k)
  
  xnew(PAN) = (xold(PAN) + dt2 * P) / (1.0 + dt2 * L)


!-> MPAN

  P = rct(72,k) * xnew(MACRO2) * xnew(NO2)
  
  L = 2.9e-11 * xnew(OH)  &
    + rct(44,k)
  
  xnew(MPAN) = (xold(MPAN) + dt2 * P) / (1.0 + dt2 * L)


!-> APINENE

  P = RCBIO(2,K)
  
  L = rct(75,k) * xnew(OH)  &
    + rct(78,k) * xnew(O3)  &
    + rct(79,k) * xnew(NO3)
  
  xnew(APINENE) = (xold(APINENE) + dt2 * P) / (1.0 + dt2 * L)


!-> TERPOOH

  P = rct(76,k) * xnew(TERPO2) * xnew(HO2)
  
  L = rct(77,k) * xnew(OH)
  
  xnew(TERPOOH) = (xold(TERPOOH) + dt2 * P) / (1.0 + dt2 * L)


!-> SO2

  P = rcemis(SO2,k)
  
  L = AQRCK(ICLOHSO2,K)*2.00e-12 * xnew(OH)  &
    + AQRCK(ICLRC1,K) * xnew(H2O2)  &
    + AQRCK(ICLRC2,K) * xnew(O3)  &
    + AQRCK(ICLRC3,K)
  
  xnew(SO2) = (xold(SO2) + dt2 * P) / (1.0 + dt2 * L)


!-> shipNOx

  P = rcemis(SHIPNOX,k)
  
  L = rct(20,k) * xnew(OH)  &
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


!-> SQT_SOA_NV

  P = rcemis(SQT_SOA_NV,k)  &
    + 0.00383* RCBIO(2,K)
  
  ! L = 0.0
  
  xnew(SQT_SOA_NV) = xold(SQT_SOA_NV) + dt2 * P


!-> POLLEN_BIRCH

  P = rcemis(POLLEN_BIRCH,k)
  
  ! L = 0.0
  
  xnew(POLLEN_BIRCH) = xold(POLLEN_BIRCH) + dt2 * P


!-> POLLEN_OLIVE

  P = rcemis(POLLEN_OLIVE,k)
  
  ! L = 0.0
  
  xnew(POLLEN_OLIVE) = xold(POLLEN_OLIVE) + dt2 * P


!-> POLLEN_ALDER

  P = rcemis(POLLEN_ALDER,k)
  
  ! L = 0.0
  
  xnew(POLLEN_ALDER) = xold(POLLEN_ALDER) + dt2 * P


!-> POLLEN_RWEED

  P = rcemis(POLLEN_RWEED,k)
  
  ! L = 0.0
  
  xnew(POLLEN_RWEED) = xold(POLLEN_RWEED) + dt2 * P


!-> POLLEN_GRASS

  P = rcemis(POLLEN_GRASS,k)
  
  ! L = 0.0
  
  xnew(POLLEN_GRASS) = xold(POLLEN_GRASS) + dt2 * P


!-> POLLEN_MUGWORT1

  P = rcemis(POLLEN_MUGWORT1,k)
  
  ! L = 0.0
  
  xnew(POLLEN_MUGWORT1) = xold(POLLEN_MUGWORT1) + dt2 * P


!-> POLLEN_MUGWORT2

  P = rcemis(POLLEN_MUGWORT2,k)
  
  ! L = 0.0
  
  xnew(POLLEN_MUGWORT2) = xold(POLLEN_MUGWORT2) + dt2 * P


!-> POLLEN_MUGWORT3

  P = rcemis(POLLEN_MUGWORT3,k)
  
  ! L = 0.0
  
  xnew(POLLEN_MUGWORT3) = xold(POLLEN_MUGWORT3) + dt2 * P


!-> POLLEN_MUGWORT4

  P = rcemis(POLLEN_MUGWORT4,k)
  
  ! L = 0.0
  
  xnew(POLLEN_MUGWORT4) = xold(POLLEN_MUGWORT4) + dt2 * P


!-> POLLEN_MUGWORT5

  P = rcemis(POLLEN_MUGWORT5,k)
  
  ! L = 0.0
  
  xnew(POLLEN_MUGWORT5) = xold(POLLEN_MUGWORT5) + dt2 * P

