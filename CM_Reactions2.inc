!-> ASOC_ng1e2

  P = rct(78,k) * xnew(OH) * xnew(ASOC_ug1)
  
  ! L = 0.0
  
  xnew(ASOC_ng1e2) = xold(ASOC_ng1e2) + dt2 * P


!-> ASOC_ug1

  P = YCOXY(0)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YCBENZ(0)* rct(47,k) * xnew(BENZENE) * xnew(OH)  &
    + YCTOL(0)* rct(48,k) * xnew(TOLUENE) * xnew(OH)  &
    + YCALK(0)* rct(39,k) * xnew(OH) * xnew(NC4H10)  &
    + YCOLE(0)* rct(73,k) * xnew(OH) * xnew(C3H6)  &
    + rct(79,k) * xnew(OH) * xnew(ASOC_ug10)
  
  L = rct(78,k) * xnew(OH)
  
  xnew(ASOC_ug1) = (xold(ASOC_ug1) + dt2 * P) / (1.0 + dt2 * L)


!-> ASOC_ug10

  P = YCOXY(1)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YCBENZ(1)* rct(47,k) * xnew(BENZENE) * xnew(OH)  &
    + YCTOL(1)* rct(48,k) * xnew(TOLUENE) * xnew(OH)  &
    + YCALK(1)* rct(39,k) * xnew(OH) * xnew(NC4H10)  &
    + YCOLE(1)* rct(73,k) * xnew(OH) * xnew(C3H6)  &
    + rct(80,k) * xnew(OH) * xnew(ASOC_ug1e2)
  
  L = rct(79,k) * xnew(OH)
  
  xnew(ASOC_ug10) = (xold(ASOC_ug10) + dt2 * P) / (1.0 + dt2 * L)


!-> ASOC_ug1e2

  P = YCOXY(2)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YCBENZ(2)* rct(47,k) * xnew(BENZENE) * xnew(OH)  &
    + YCTOL(2)* rct(48,k) * xnew(TOLUENE) * xnew(OH)  &
    + YCALK(2)* rct(39,k) * xnew(OH) * xnew(NC4H10)  &
    + YCOLE(2)* rct(73,k) * xnew(OH) * xnew(C3H6)  &
    + rct(81,k) * xnew(OH) * xnew(ASOC_ug1e3)
  
  L = rct(80,k) * xnew(OH)
  
  xnew(ASOC_ug1e2) = (xold(ASOC_ug1e2) + dt2 * P) / (1.0 + dt2 * L)


!-> ASOC_ug1e3

  P = YCOXY(3)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YCBENZ(3)* rct(47,k) * xnew(BENZENE) * xnew(OH)  &
    + YCTOL(3)* rct(48,k) * xnew(TOLUENE) * xnew(OH)  &
    + YCALK(3)* rct(39,k) * xnew(OH) * xnew(NC4H10)  &
    + YCOLE(3)* rct(73,k) * xnew(OH) * xnew(C3H6)
  
  L = rct(81,k) * xnew(OH)
  
  xnew(ASOC_ug1e3) = (xold(ASOC_ug1e3) + dt2 * P) / (1.0 + dt2 * L)


!-> non_C_ASOA_ng1e2

  P = 0.9* rct(78,k) * xnew(OH) * xnew(ASOC_ug1)  &
    + 1.075* rct(82,k) * xnew(OH) * xnew(non_C_ASOA_ug1)
  
  ! L = 0.0
  
  xnew(non_C_ASOA_ng1e2) = xold(non_C_ASOA_ng1e2) + dt2 * P


!-> non_C_ASOA_ug1

  P = YNOXY(0)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YNBENZ(0)* rct(47,k) * xnew(BENZENE) * xnew(OH)  &
    + YNTOL(0)* rct(48,k) * xnew(TOLUENE) * xnew(OH)  &
    + YNALK(0)* rct(39,k) * xnew(OH) * xnew(NC4H10)  &
    + YNOLE(0)* rct(73,k) * xnew(OH) * xnew(C3H6)  &
    + 0.9* rct(79,k) * xnew(OH) * xnew(ASOC_ug10)  &
    + 1.075* rct(83,k) * xnew(OH) * xnew(non_C_ASOA_ug10)
  
  L = rct(82,k) * xnew(OH)
  
  xnew(non_C_ASOA_ug1) = (xold(non_C_ASOA_ug1) + dt2 * P) / (1.0 + dt2 * L)


!-> non_C_ASOA_ug10

  P = YNOXY(1)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YNBENZ(1)* rct(47,k) * xnew(BENZENE) * xnew(OH)  &
    + YNTOL(1)* rct(48,k) * xnew(TOLUENE) * xnew(OH)  &
    + YNALK(1)* rct(39,k) * xnew(OH) * xnew(NC4H10)  &
    + YNOLE(1)* rct(73,k) * xnew(OH) * xnew(C3H6)  &
    + 0.9* rct(80,k) * xnew(OH) * xnew(ASOC_ug1e2)  &
    + 1.075* rct(84,k) * xnew(OH) * xnew(non_C_ASOA_ug1e2)
  
  L = rct(83,k) * xnew(OH)
  
  xnew(non_C_ASOA_ug10) = (xold(non_C_ASOA_ug10) + dt2 * P) / (1.0 + dt2 * L)


!-> non_C_ASOA_ug1e2

  P = YNOXY(2)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YNBENZ(2)* rct(47,k) * xnew(BENZENE) * xnew(OH)  &
    + YNTOL(2)* rct(48,k) * xnew(TOLUENE) * xnew(OH)  &
    + YNALK(2)* rct(39,k) * xnew(OH) * xnew(NC4H10)  &
    + YNOLE(2)* rct(73,k) * xnew(OH) * xnew(C3H6)  &
    + 0.9* rct(81,k) * xnew(OH) * xnew(ASOC_ug1e3)  &
    + 1.075* rct(85,k) * xnew(OH) * xnew(non_C_ASOA_ug1e3)
  
  L = rct(84,k) * xnew(OH)
  
  xnew(non_C_ASOA_ug1e2) = (xold(non_C_ASOA_ug1e2) + dt2 * P) / (1.0 + dt2 * L)


!-> non_C_ASOA_ug1e3

  P = YNOXY(3)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YNBENZ(3)* rct(47,k) * xnew(BENZENE) * xnew(OH)  &
    + YNTOL(3)* rct(48,k) * xnew(TOLUENE) * xnew(OH)  &
    + YNALK(3)* rct(39,k) * xnew(OH) * xnew(NC4H10)  &
    + YNOLE(3)* rct(73,k) * xnew(OH) * xnew(C3H6)
  
  L = rct(85,k) * xnew(OH)
  
  xnew(non_C_ASOA_ug1e3) = (xold(non_C_ASOA_ug1e3) + dt2 * P) / (1.0 + dt2 * L)


!-> BSOC_ng1e2

  P = rct(86,k) * xnew(OH) * xnew(BSOC_ug1)
  
  ! L = 0.0
  
  xnew(BSOC_ng1e2) = xold(BSOC_ng1e2) + dt2 * P


!-> BSOC_ug1

  P = YCISOP(0)* rct(55,k) * xnew(OH) * xnew(C5H8)  &
    + YCTERP(0)* rct(75,k) * xnew(OH) * xnew(APINENE)  &
    + YCTERP(0)* rct(76,k) * xnew(O3) * xnew(APINENE)  &
    + YCTERP(0)* rct(77,k) * xnew(NO3) * xnew(APINENE)  &
    + rct(87,k) * xnew(OH) * xnew(BSOC_ug10)
  
  L = rct(86,k) * xnew(OH)
  
  xnew(BSOC_ug1) = (xold(BSOC_ug1) + dt2 * P) / (1.0 + dt2 * L)


!-> BSOC_ug10

  P = YCISOP(1)* rct(55,k) * xnew(OH) * xnew(C5H8)  &
    + YCTERP(1)* rct(75,k) * xnew(OH) * xnew(APINENE)  &
    + YCTERP(1)* rct(76,k) * xnew(O3) * xnew(APINENE)  &
    + YCTERP(1)* rct(77,k) * xnew(NO3) * xnew(APINENE)  &
    + rct(88,k) * xnew(OH) * xnew(BSOC_ug1e2)
  
  L = rct(87,k) * xnew(OH)
  
  xnew(BSOC_ug10) = (xold(BSOC_ug10) + dt2 * P) / (1.0 + dt2 * L)


!-> BSOC_ug1e2

  P = YCISOP(2)* rct(55,k) * xnew(OH) * xnew(C5H8)  &
    + YCTERP(2)* rct(75,k) * xnew(OH) * xnew(APINENE)  &
    + YCTERP(2)* rct(76,k) * xnew(O3) * xnew(APINENE)  &
    + YCTERP(2)* rct(77,k) * xnew(NO3) * xnew(APINENE)  &
    + rct(89,k) * xnew(OH) * xnew(BSOC_ug1e3)
  
  L = rct(88,k) * xnew(OH)
  
  xnew(BSOC_ug1e2) = (xold(BSOC_ug1e2) + dt2 * P) / (1.0 + dt2 * L)


!-> BSOC_ug1e3

  P = YCISOP(3)* rct(55,k) * xnew(OH) * xnew(C5H8)  &
    + YCTERP(3)* rct(75,k) * xnew(OH) * xnew(APINENE)  &
    + YCTERP(3)* rct(76,k) * xnew(O3) * xnew(APINENE)  &
    + YCTERP(3)* rct(77,k) * xnew(NO3) * xnew(APINENE)
  
  L = rct(89,k) * xnew(OH)
  
  xnew(BSOC_ug1e3) = (xold(BSOC_ug1e3) + dt2 * P) / (1.0 + dt2 * L)


!-> non_C_BSOA_ng1e2

  P = 0.9* rct(86,k) * xnew(OH) * xnew(BSOC_ug1)  &
    + 1.075* rct(90,k) * xnew(OH) * xnew(non_C_BSOA_ug1)
  
  ! L = 0.0
  
  xnew(non_C_BSOA_ng1e2) = xold(non_C_BSOA_ng1e2) + dt2 * P


!-> non_C_BSOA_ug1

  P = YNISOP(0)* rct(55,k) * xnew(OH) * xnew(C5H8)  &
    + YNTERP(0)* rct(75,k) * xnew(OH) * xnew(APINENE)  &
    + YNTERP(0)* rct(76,k) * xnew(O3) * xnew(APINENE)  &
    + YNTERP(0)* rct(77,k) * xnew(NO3) * xnew(APINENE)  &
    + 0.9* rct(87,k) * xnew(OH) * xnew(BSOC_ug10)  &
    + 1.075* rct(91,k) * xnew(OH) * xnew(non_C_BSOA_ug10)
  
  L = rct(90,k) * xnew(OH)
  
  xnew(non_C_BSOA_ug1) = (xold(non_C_BSOA_ug1) + dt2 * P) / (1.0 + dt2 * L)


!-> non_C_BSOA_ug10

  P = YNISOP(1)* rct(55,k) * xnew(OH) * xnew(C5H8)  &
    + YNTERP(1)* rct(75,k) * xnew(OH) * xnew(APINENE)  &
    + YNTERP(1)* rct(76,k) * xnew(O3) * xnew(APINENE)  &
    + YNTERP(1)* rct(77,k) * xnew(NO3) * xnew(APINENE)  &
    + 0.9* rct(88,k) * xnew(OH) * xnew(BSOC_ug1e2)  &
    + 1.075* rct(92,k) * xnew(OH) * xnew(non_C_BSOA_ug1e2)
  
  L = rct(91,k) * xnew(OH)
  
  xnew(non_C_BSOA_ug10) = (xold(non_C_BSOA_ug10) + dt2 * P) / (1.0 + dt2 * L)


!-> non_C_BSOA_ug1e2

  P = YNISOP(2)* rct(55,k) * xnew(OH) * xnew(C5H8)  &
    + YNTERP(2)* rct(75,k) * xnew(OH) * xnew(APINENE)  &
    + YNTERP(2)* rct(76,k) * xnew(O3) * xnew(APINENE)  &
    + YNTERP(2)* rct(77,k) * xnew(NO3) * xnew(APINENE)  &
    + 0.9* rct(89,k) * xnew(OH) * xnew(BSOC_ug1e3)  &
    + 1.075* rct(93,k) * xnew(OH) * xnew(non_C_BSOA_ug1e3)
  
  L = rct(92,k) * xnew(OH)
  
  xnew(non_C_BSOA_ug1e2) = (xold(non_C_BSOA_ug1e2) + dt2 * P) / (1.0 + dt2 * L)


!-> non_C_BSOA_ug1e3

  P = YNISOP(3)* rct(55,k) * xnew(OH) * xnew(C5H8)  &
    + YNTERP(3)* rct(75,k) * xnew(OH) * xnew(APINENE)  &
    + YNTERP(3)* rct(76,k) * xnew(O3) * xnew(APINENE)  &
    + YNTERP(3)* rct(77,k) * xnew(NO3) * xnew(APINENE)
  
  L = rct(93,k) * xnew(OH)
  
  xnew(non_C_BSOA_ug1e3) = (xold(non_C_BSOA_ug1e3) + dt2 * P) / (1.0 + dt2 * L)


!-> ffuel_ng10

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for ffuel_ng10! xnew(ffuel_ng10) = max(0.0, xold(ffuel_ng10))


!-> woodOA_ng10

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for woodOA_ng10! xnew(woodOA_ng10) = max(0.0, xold(woodOA_ng10))


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

  P = rct(102,k) * xnew(HNO3)  &
    + rct(103,k) * xnew(HNO3)
  
  ! L = 0.0
  
  xnew(NO3_c) = xold(NO3_c) + dt2 * P


!-> NH4_f

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for NH4_f! xnew(NH4_f) = max(0.0, xold(NH4_f))


!-> POM_f_wood

  P = rcemis(POM_F_WOOD,k)
  
  ! L = 0.0
  
  xnew(POM_f_wood) = xold(POM_f_wood) + dt2 * P


!-> POM_c_wood

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for POM_c_wood! xnew(POM_c_wood) = max(0.0, xold(POM_c_wood))


!-> POM_f_ffuel

  P = rcemis(POM_F_FFUEL,k)
  
  ! L = 0.0
  
  xnew(POM_f_ffuel) = xold(POM_f_ffuel) + dt2 * P


!-> POM_c_ffuel

  P = rcemis(POM_C_FFUEL,k)
  
  ! L = 0.0
  
  xnew(POM_c_ffuel) = xold(POM_c_ffuel) + dt2 * P


!-> EC_f_wood_new

  P = rcemis(EC_F_WOOD_new,k)
  
  L = rct(94,k)
  
  xnew(EC_f_wood_new) = (xold(EC_f_wood_new) + dt2 * P) / (1.0 + dt2 * L)


!-> EC_f_wood_age

  P = rcemis(EC_F_WOOD_age,k)  &
    + rct(94,k) * xnew(EC_f_wood_new)
  
  ! L = 0.0
  
  xnew(EC_f_wood_age) = xold(EC_f_wood_age) + dt2 * P


!-> EC_c_wood

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for EC_c_wood! xnew(EC_c_wood) = max(0.0, xold(EC_c_wood))


!-> EC_f_ffuel_new

  P = rcemis(EC_F_FFUEL_new,k)
  
  L = rct(94,k)
  
  xnew(EC_f_ffuel_new) = (xold(EC_f_ffuel_new) + dt2 * P) / (1.0 + dt2 * L)


!-> EC_f_ffuel_age

  P = rcemis(EC_F_FFUEL_age,k)  &
    + rct(94,k) * xnew(EC_f_ffuel_new)
  
  ! L = 0.0
  
  xnew(EC_f_ffuel_age) = xold(EC_f_ffuel_age) + dt2 * P


!-> EC_c_ffuel

  P = rcemis(EC_C_FFUEL,k)
  
  ! L = 0.0
  
  xnew(EC_c_ffuel) = xold(EC_c_ffuel) + dt2 * P


!-> pSO4f

  P = rcemis(PSO4F,k)
  
  ! L = 0.0
  
  xnew(pSO4f) = xold(pSO4f) + dt2 * P


!-> pSO4c

  P = rcemis(PSO4C,k)
  
  ! L = 0.0
  
  xnew(pSO4c) = xold(pSO4c) + dt2 * P


!-> remPPM25

  P = rcemis(REMPPM25,k)
  
  ! L = 0.0
  
  xnew(remPPM25) = xold(remPPM25) + dt2 * P


!-> remPPM_c

  P = rcemis(REMPPM_C,k)
  
  ! L = 0.0
  
  xnew(remPPM_c) = xold(remPPM_c) + dt2 * P


!-> OM25_bgnd

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for OM25_bgnd! xnew(OM25_bgnd) = max(0.0, xold(OM25_bgnd))


!-> OM25_p

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for OM25_p! xnew(OM25_p) = max(0.0, xold(OM25_p))


!-> Ash_f

  P = rcemis(ASH_F,K)
  
  ! L = 0.0
  
  xnew(Ash_f) = xold(Ash_f) + dt2 * P


!-> Ash_c

  P = rcemis(ASH_C,K)
  
  ! L = 0.0
  
  xnew(Ash_c) = xold(Ash_c) + dt2 * P


!-> ffire_CO

  P = rcemis(FFIRE_CO,k)
  
  L = rct(105,k) * xnew(OH)
  
  xnew(ffire_CO) = (xold(ffire_CO) + dt2 * P) / (1.0 + dt2 * L)


!-> ffire_OM

  P = rcemis(FFIRE_OM,k)
  
  ! L = 0.0
  
  xnew(ffire_OM) = xold(ffire_OM) + dt2 * P


!-> ffire_BC

  P = rcemis(FFIRE_BC,k)
  
  ! L = 0.0
  
  xnew(ffire_BC) = xold(ffire_BC) + dt2 * P


!-> ffire_remPPM25

  P = rcemis(FFIRE_REMPPM25,k)
  
  ! L = 0.0
  
  xnew(ffire_remPPM25) = xold(ffire_remPPM25) + dt2 * P


!-> SeaSalt_f

  P = rcemis(SEASALT_F,K)
  
  ! L = 0.0
  
  xnew(SeaSalt_f) = xold(SeaSalt_f) + dt2 * P


!-> SeaSalt_c

  P = rcemis(SEASALT_C,K)
  
  ! L = 0.0
  
  xnew(SeaSalt_c) = xold(SeaSalt_c) + dt2 * P


