!-> ASOC_ng100

  P = rct(100,k) * xnew(OH) * xnew(ASOC_ug1)
  
  ! L = 0.0
  
  xnew(ASOC_ng100) = xold(ASOC_ng100) + dt2 * P


!-> ASOC_ug1

  P = YCOXY(0)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YCALK(0)* rct(44,k) * xnew(OH) * xnew(NC4H10)  &
    + YCOLE(0)* rct(82,k) * xnew(OH) * xnew(C3H6)  &
    + rct(101,k) * xnew(OH) * xnew(ASOC_ug10)
  
  L = rct(100,k) * xnew(OH)
  
  xnew(ASOC_ug1) = (xold(ASOC_ug1) + dt2 * P) / (1.0 + dt2 * L)


!-> ASOC_ug10

  P = YCOXY(1)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YCALK(1)* rct(44,k) * xnew(OH) * xnew(NC4H10)  &
    + YCOLE(1)* rct(82,k) * xnew(OH) * xnew(C3H6)  &
    + rct(102,k) * xnew(OH) * xnew(ASOC_ug1e2)
  
  L = rct(101,k) * xnew(OH)
  
  xnew(ASOC_ug10) = (xold(ASOC_ug10) + dt2 * P) / (1.0 + dt2 * L)


!-> ASOC_ug1e2

  P = YCOXY(2)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YCALK(2)* rct(44,k) * xnew(OH) * xnew(NC4H10)  &
    + YCOLE(2)* rct(82,k) * xnew(OH) * xnew(C3H6)  &
    + rct(103,k) * xnew(OH) * xnew(ASOC_ug1e3)
  
  L = rct(102,k) * xnew(OH)
  
  xnew(ASOC_ug1e2) = (xold(ASOC_ug1e2) + dt2 * P) / (1.0 + dt2 * L)


!-> ASOC_ug1e3

  P = YCOXY(3)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YCALK(3)* rct(44,k) * xnew(OH) * xnew(NC4H10)  &
    + YCOLE(3)* rct(82,k) * xnew(OH) * xnew(C3H6)
  
  L = rct(103,k) * xnew(OH)
  
  xnew(ASOC_ug1e3) = (xold(ASOC_ug1e3) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_ASOA_ng100

  P = 0.9* rct(100,k) * xnew(OH) * xnew(ASOC_ug1)  &
    + 1.075* rct(104,k) * xnew(OH) * xnew(NON_C_ASOA_ug1)
  
  ! L = 0.0
  
  xnew(NON_C_ASOA_ng100) = xold(NON_C_ASOA_ng100) + dt2 * P


!-> NON_C_ASOA_ug1

  P = YNOXY(0)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YNALK(0)* rct(44,k) * xnew(OH) * xnew(NC4H10)  &
    + YNOLE(0)* rct(82,k) * xnew(OH) * xnew(C3H6)  &
    + 0.9* rct(101,k) * xnew(OH) * xnew(ASOC_ug10)  &
    + 1.075* rct(105,k) * xnew(OH) * xnew(NON_C_ASOA_ug10)
  
  L = rct(104,k) * xnew(OH)
  
  xnew(NON_C_ASOA_ug1) = (xold(NON_C_ASOA_ug1) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_ASOA_ug10

  P = YNOXY(1)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YNALK(1)* rct(44,k) * xnew(OH) * xnew(NC4H10)  &
    + YNOLE(1)* rct(82,k) * xnew(OH) * xnew(C3H6)  &
    + 0.9* rct(102,k) * xnew(OH) * xnew(ASOC_ug1e2)  &
    + 1.075* rct(106,k) * xnew(OH) * xnew(NON_C_ASOA_ug1e2)
  
  L = rct(105,k) * xnew(OH)
  
  xnew(NON_C_ASOA_ug10) = (xold(NON_C_ASOA_ug10) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_ASOA_ug1e2

  P = YNOXY(2)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YNALK(2)* rct(44,k) * xnew(OH) * xnew(NC4H10)  &
    + YNOLE(2)* rct(82,k) * xnew(OH) * xnew(C3H6)  &
    + 0.9* rct(103,k) * xnew(OH) * xnew(ASOC_ug1e3)  &
    + 1.075* rct(107,k) * xnew(OH) * xnew(NON_C_ASOA_ug1e3)
  
  L = rct(106,k) * xnew(OH)
  
  xnew(NON_C_ASOA_ug1e2) = (xold(NON_C_ASOA_ug1e2) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_ASOA_ug1e3

  P = YNOXY(3)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YNALK(3)* rct(44,k) * xnew(OH) * xnew(NC4H10)  &
    + YNOLE(3)* rct(82,k) * xnew(OH) * xnew(C3H6)
  
  L = rct(107,k) * xnew(OH)
  
  xnew(NON_C_ASOA_ug1e3) = (xold(NON_C_ASOA_ug1e3) + dt2 * P) / (1.0 + dt2 * L)


!-> BSOC_ng100

  P = rct(108,k) * xnew(OH) * xnew(BSOC_ug1)
  
  ! L = 0.0
  
  xnew(BSOC_ng100) = xold(BSOC_ng100) + dt2 * P


!-> BSOC_ug1

  P = YCISOP(0)* rct(55,k) * xnew(OH) * xnew(C5H8)  &
    + YCTERP(0)* rct(84,k) * xnew(OH) * xnew(APINENE)  &
    + YCTERP(0)* rct(85,k) * xnew(OH) * xnew(BPINENE)  &
    + YCTERP(0)* rct(86,k) * xnew(OH) * xnew(XTERP)  &
    + YCTERP(0)* rct(89,k) * xnew(O3) * xnew(APINENE)  &
    + YCTERP(0)* rct(90,k) * xnew(O3) * xnew(BPINENE)  &
    + YCTERP(0)* rct(91,k) * xnew(O3) * xnew(XTERP)  &
    + YCTERP(0)* rct(99,k) * xnew(NO3) * xnew(APINENE)  &
    + YCTERP(0)* 2.51e-12 * xnew(NO3) * xnew(BPINENE)  &
    + YCTERP(0)* rct(93,k) * xnew(NO3) * xnew(XTERP)  &
    + rct(109,k) * xnew(OH) * xnew(BSOC_ug10)
  
  L = rct(108,k) * xnew(OH)
  
  xnew(BSOC_ug1) = (xold(BSOC_ug1) + dt2 * P) / (1.0 + dt2 * L)


!-> BSOC_ug10

  P = YCISOP(1)* rct(55,k) * xnew(OH) * xnew(C5H8)  &
    + YCTERP(1)* rct(84,k) * xnew(OH) * xnew(APINENE)  &
    + YCTERP(1)* rct(85,k) * xnew(OH) * xnew(BPINENE)  &
    + YCTERP(1)* rct(86,k) * xnew(OH) * xnew(XTERP)  &
    + YCTERP(1)* rct(89,k) * xnew(O3) * xnew(APINENE)  &
    + YCTERP(1)* rct(90,k) * xnew(O3) * xnew(BPINENE)  &
    + YCTERP(1)* rct(91,k) * xnew(O3) * xnew(XTERP)  &
    + YCTERP(1)* rct(99,k) * xnew(NO3) * xnew(APINENE)  &
    + YCTERP(1)* 2.51e-12 * xnew(NO3) * xnew(BPINENE)  &
    + YCTERP(1)* rct(93,k) * xnew(NO3) * xnew(XTERP)  &
    + rct(110,k) * xnew(OH) * xnew(BSOC_ug1e2)
  
  L = rct(109,k) * xnew(OH)
  
  xnew(BSOC_ug10) = (xold(BSOC_ug10) + dt2 * P) / (1.0 + dt2 * L)


!-> BSOC_ug1e2

  P = YCISOP(2)* rct(55,k) * xnew(OH) * xnew(C5H8)  &
    + YCTERP(2)* rct(84,k) * xnew(OH) * xnew(APINENE)  &
    + YCTERP(2)* rct(85,k) * xnew(OH) * xnew(BPINENE)  &
    + YCTERP(2)* rct(86,k) * xnew(OH) * xnew(XTERP)  &
    + YCTERP(2)* rct(89,k) * xnew(O3) * xnew(APINENE)  &
    + YCTERP(2)* rct(90,k) * xnew(O3) * xnew(BPINENE)  &
    + YCTERP(2)* rct(91,k) * xnew(O3) * xnew(XTERP)  &
    + YCTERP(2)* rct(99,k) * xnew(NO3) * xnew(APINENE)  &
    + YCTERP(2)* 2.51e-12 * xnew(NO3) * xnew(BPINENE)  &
    + YCTERP(2)* rct(93,k) * xnew(NO3) * xnew(XTERP)  &
    + rct(111,k) * xnew(OH) * xnew(BSOC_ug1e3)
  
  L = rct(110,k) * xnew(OH)
  
  xnew(BSOC_ug1e2) = (xold(BSOC_ug1e2) + dt2 * P) / (1.0 + dt2 * L)


!-> BSOC_ug1e3

  P = YCISOP(3)* rct(55,k) * xnew(OH) * xnew(C5H8)  &
    + YCTERP(3)* rct(84,k) * xnew(OH) * xnew(APINENE)  &
    + YCTERP(3)* rct(85,k) * xnew(OH) * xnew(BPINENE)  &
    + YCTERP(3)* rct(86,k) * xnew(OH) * xnew(XTERP)  &
    + YCTERP(3)* rct(89,k) * xnew(O3) * xnew(APINENE)  &
    + YCTERP(3)* rct(90,k) * xnew(O3) * xnew(BPINENE)  &
    + YCTERP(3)* rct(91,k) * xnew(O3) * xnew(XTERP)  &
    + YCTERP(3)* rct(99,k) * xnew(NO3) * xnew(APINENE)  &
    + YCTERP(3)* 2.51e-12 * xnew(NO3) * xnew(BPINENE)  &
    + YCTERP(3)* rct(93,k) * xnew(NO3) * xnew(XTERP)
  
  L = rct(111,k) * xnew(OH)
  
  xnew(BSOC_ug1e3) = (xold(BSOC_ug1e3) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_BSOA_ng100

  P = 0.9* rct(108,k) * xnew(OH) * xnew(BSOC_ug1)  &
    + 1.075* rct(112,k) * xnew(OH) * xnew(NON_C_BSOA_ug1)
  
  ! L = 0.0
  
  xnew(NON_C_BSOA_ng100) = xold(NON_C_BSOA_ng100) + dt2 * P


!-> NON_C_BSOA_ug1

  P = YNISOP(0)* rct(55,k) * xnew(OH) * xnew(C5H8)  &
    + YNTERP(0)* rct(84,k) * xnew(OH) * xnew(APINENE)  &
    + YNTERP(0)* rct(85,k) * xnew(OH) * xnew(BPINENE)  &
    + YNTERP(0)* rct(86,k) * xnew(OH) * xnew(XTERP)  &
    + YNTERP(0)* rct(89,k) * xnew(O3) * xnew(APINENE)  &
    + YNTERP(0)* rct(90,k) * xnew(O3) * xnew(BPINENE)  &
    + YNTERP(0)* rct(91,k) * xnew(O3) * xnew(XTERP)  &
    + YNTERP(0)* rct(99,k) * xnew(NO3) * xnew(APINENE)  &
    + YNTERP(0)* 2.51e-12 * xnew(NO3) * xnew(BPINENE)  &
    + YNTERP(0)* rct(93,k) * xnew(NO3) * xnew(XTERP)  &
    + 0.9* rct(109,k) * xnew(OH) * xnew(BSOC_ug10)  &
    + 1.075* rct(113,k) * xnew(OH) * xnew(NON_C_BSOA_ug10)
  
  L = rct(112,k) * xnew(OH)
  
  xnew(NON_C_BSOA_ug1) = (xold(NON_C_BSOA_ug1) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_BSOA_ug10

  P = YNISOP(1)* rct(55,k) * xnew(OH) * xnew(C5H8)  &
    + YNTERP(1)* rct(84,k) * xnew(OH) * xnew(APINENE)  &
    + YNTERP(1)* rct(85,k) * xnew(OH) * xnew(BPINENE)  &
    + YNTERP(1)* rct(86,k) * xnew(OH) * xnew(XTERP)  &
    + YNTERP(1)* rct(89,k) * xnew(O3) * xnew(APINENE)  &
    + YNTERP(1)* rct(90,k) * xnew(O3) * xnew(BPINENE)  &
    + YNTERP(1)* rct(91,k) * xnew(O3) * xnew(XTERP)  &
    + YNTERP(1)* rct(99,k) * xnew(NO3) * xnew(APINENE)  &
    + YNTERP(1)* 2.51e-12 * xnew(NO3) * xnew(BPINENE)  &
    + YNTERP(1)* rct(93,k) * xnew(NO3) * xnew(XTERP)  &
    + 0.9* rct(110,k) * xnew(OH) * xnew(BSOC_ug1e2)  &
    + 1.075* rct(114,k) * xnew(OH) * xnew(NON_C_BSOA_ug1e2)
  
  L = rct(113,k) * xnew(OH)
  
  xnew(NON_C_BSOA_ug10) = (xold(NON_C_BSOA_ug10) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_BSOA_ug1e2

  P = YNISOP(2)* rct(55,k) * xnew(OH) * xnew(C5H8)  &
    + YNTERP(2)* rct(84,k) * xnew(OH) * xnew(APINENE)  &
    + YNTERP(2)* rct(85,k) * xnew(OH) * xnew(BPINENE)  &
    + YNTERP(2)* rct(86,k) * xnew(OH) * xnew(XTERP)  &
    + YNTERP(2)* rct(89,k) * xnew(O3) * xnew(APINENE)  &
    + YNTERP(2)* rct(90,k) * xnew(O3) * xnew(BPINENE)  &
    + YNTERP(2)* rct(91,k) * xnew(O3) * xnew(XTERP)  &
    + YNTERP(2)* rct(99,k) * xnew(NO3) * xnew(APINENE)  &
    + YNTERP(2)* 2.51e-12 * xnew(NO3) * xnew(BPINENE)  &
    + YNTERP(2)* rct(93,k) * xnew(NO3) * xnew(XTERP)  &
    + 0.9* rct(111,k) * xnew(OH) * xnew(BSOC_ug1e3)  &
    + 1.075* rct(115,k) * xnew(OH) * xnew(NON_C_BSOA_ug1e3)
  
  L = rct(114,k) * xnew(OH)
  
  xnew(NON_C_BSOA_ug1e2) = (xold(NON_C_BSOA_ug1e2) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_BSOA_ug1e3

  P = YNISOP(3)* rct(55,k) * xnew(OH) * xnew(C5H8)  &
    + YNTERP(3)* rct(84,k) * xnew(OH) * xnew(APINENE)  &
    + YNTERP(3)* rct(85,k) * xnew(OH) * xnew(BPINENE)  &
    + YNTERP(3)* rct(86,k) * xnew(OH) * xnew(XTERP)  &
    + YNTERP(3)* rct(89,k) * xnew(O3) * xnew(APINENE)  &
    + YNTERP(3)* rct(90,k) * xnew(O3) * xnew(BPINENE)  &
    + YNTERP(3)* rct(91,k) * xnew(O3) * xnew(XTERP)  &
    + YNTERP(3)* rct(99,k) * xnew(NO3) * xnew(APINENE)  &
    + YNTERP(3)* 2.51e-12 * xnew(NO3) * xnew(BPINENE)  &
    + YNTERP(3)* rct(93,k) * xnew(NO3) * xnew(XTERP)
  
  L = rct(115,k) * xnew(OH)
  
  xnew(NON_C_BSOA_ug1e3) = (xold(NON_C_BSOA_ug1e3) + dt2 * P) / (1.0 + dt2 * L)


!-> FFFUEL_ng10

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for FFFUEL_ng10! xnew(FFFUEL_ng10) = max(0.0, xold(FFFUEL_ng10))


!-> WOODOA_ng10

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for WOODOA_ng10! xnew(WOODOA_ng10) = max(0.0, xold(WOODOA_ng10))


!-> SO4

  P = rcemis(SO4,k)  &
    + AQRCK(ICLOHSO2,K)*2.00e-12 * xnew(OH) * xnew(SO2)  &
    + AQRCK(ICLRC1,K) * xnew(SO2) * xnew(H2O2)  &
    + AQRCK(ICLRC2,K) * xnew(SO2) * xnew(O3)  &
    + AQRCK(ICLRC3,K) * xnew(SO2)
  
  ! L = 0.0
  
  xnew(SO4) = xold(SO4) + dt2 * P


!-> NH3

  P = rcemis(NH3,k)
  
  ! L = 0.0
  
  xnew(NH3) = xold(NH3) + dt2 * P


!-> NO3_f

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for NO3_f! xnew(NO3_f) = max(0.0, xold(NO3_f))


!-> NO3_c

  P = rct(96,k) * xnew(HNO3)  &
    + rct(97,k) * xnew(HNO3)
  
  ! L = 0.0
  
  xnew(NO3_c) = xold(NO3_c) + dt2 * P


!-> NH4_f

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for NH4_f! xnew(NH4_f) = max(0.0, xold(NH4_f))


!-> POM_f_WOOD

  P = rcemis(POM_F_WOOD,k)
  
  ! L = 0.0
  
  xnew(POM_f_WOOD) = xold(POM_f_WOOD) + dt2 * P


!-> POM_f_FFUEL

  P = rcemis(POM_F_FFUEL,k)
  
  ! L = 0.0
  
  xnew(POM_f_FFUEL) = xold(POM_f_FFUEL) + dt2 * P


!-> POM_c_FFUEL

  P = rcemis(POM_C_FFUEL,k)
  
  ! L = 0.0
  
  xnew(POM_c_FFUEL) = xold(POM_c_FFUEL) + dt2 * P


!-> EC_f_WOOD_new

  P = rcemis(EC_F_WOOD_NEW,k)
  
  L = rct(116,k)
  
  xnew(EC_f_WOOD_new) = (xold(EC_f_WOOD_new) + dt2 * P) / (1.0 + dt2 * L)


!-> EC_f_WOOD_age

  P = rcemis(EC_F_WOOD_AGE,k)  &
    + rct(116,k) * xnew(EC_f_WOOD_new)
  
  ! L = 0.0
  
  xnew(EC_f_WOOD_age) = xold(EC_f_WOOD_age) + dt2 * P


!-> EC_c_WOOD

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for EC_c_WOOD! xnew(EC_c_WOOD) = max(0.0, xold(EC_c_WOOD))


!-> EC_f_FFUEL_new

  P = rcemis(EC_F_FFUEL_NEW,k)
  
  L = rct(116,k)
  
  xnew(EC_f_FFUEL_new) = (xold(EC_f_FFUEL_new) + dt2 * P) / (1.0 + dt2 * L)


!-> EC_f_FFUEL_age

  P = rcemis(EC_F_FFUEL_AGE,k)  &
    + rct(116,k) * xnew(EC_f_FFUEL_new)
  
  ! L = 0.0
  
  xnew(EC_f_FFUEL_age) = xold(EC_f_FFUEL_age) + dt2 * P


!-> EC_c_FFUEL

  P = rcemis(EC_C_FFUEL,k)
  
  ! L = 0.0
  
  xnew(EC_c_FFUEL) = xold(EC_c_FFUEL) + dt2 * P


!-> REMPPM25

  P = rcemis(REMPPM25,k)
  
  ! L = 0.0
  
  xnew(REMPPM25) = xold(REMPPM25) + dt2 * P


!-> REMPPM_c

  P = rcemis(REMPPM_C,k)
  
  ! L = 0.0
  
  xnew(REMPPM_c) = xold(REMPPM_c) + dt2 * P


!-> OM25_BGND

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for OM25_BGND! xnew(OM25_BGND) = max(0.0, xold(OM25_BGND))


!-> OM25_p

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for OM25_p! xnew(OM25_p) = max(0.0, xold(OM25_p))


!-> FFIRE_CO

  P = rcemis(FFIRE_CO,k)
  
  L = rct(117,k) * xnew(OH)
  
  xnew(FFIRE_CO) = (xold(FFIRE_CO) + dt2 * P) / (1.0 + dt2 * L)


!-> FFIRE_OM

  P = rcemis(FFIRE_OM,k)
  
  ! L = 0.0
  
  xnew(FFIRE_OM) = xold(FFIRE_OM) + dt2 * P


!-> FFIRE_BC

  P = rcemis(FFIRE_BC,k)
  
  ! L = 0.0
  
  xnew(FFIRE_BC) = xold(FFIRE_BC) + dt2 * P


!-> FFIRE_REMPPM25

  P = rcemis(FFIRE_REMPPM25,k)
  
  ! L = 0.0
  
  xnew(FFIRE_REMPPM25) = xold(FFIRE_REMPPM25) + dt2 * P


!-> SeaSalt_f

  P = rcemis(SEASALT_F,K)
  
  ! L = 0.0
  
  xnew(SeaSalt_f) = xold(SeaSalt_f) + dt2 * P


!-> SeaSalt_c

  P = rcemis(SEASALT_C,K)
  
  ! L = 0.0
  
  xnew(SeaSalt_c) = xold(SeaSalt_c) + dt2 * P


!-> Ash_f

  P = rcemis(ASH_F,K)
  
  ! L = 0.0
  
  xnew(Ash_f) = xold(Ash_f) + dt2 * P


!-> Ash_c

  P = rcemis(ASH_C,K)
  
  ! L = 0.0
  
  xnew(Ash_c) = xold(Ash_c) + dt2 * P


