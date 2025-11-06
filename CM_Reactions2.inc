
!-> ASOC_ug1e3

  P = YCOXY(3)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YCBENZ(3)* rct(56,k) * xnew(BENZENE) * xnew(OH)  &
    + YCTOL(3)* rct(57,k) * xnew(TOLUENE) * xnew(OH)  &
    + YCALK(3)* rct(46,k) * xnew(OH) * xnew(NC4H10)  &
    + YCOLE(3)* rct(52,k) * xnew(OH) * xnew(C3H6)
  
  L = rct(83,k) * xnew(OH)
  
  xnew(ASOC_ug1e3) = (xold(ASOC_ug1e3) + dt2 * P) / (1.0 + dt2 * L)


!-> ASOC_ug1e2

  P = YCOXY(2)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YCBENZ(2)* rct(56,k) * xnew(BENZENE) * xnew(OH)  &
    + YCTOL(2)* rct(57,k) * xnew(TOLUENE) * xnew(OH)  &
    + YCALK(2)* rct(46,k) * xnew(OH) * xnew(NC4H10)  &
    + YCOLE(2)* rct(52,k) * xnew(OH) * xnew(C3H6)  &
    + rct(83,k) * xnew(OH) * xnew(ASOC_ug1e3)
  
  L = rct(82,k) * xnew(OH)
  
  xnew(ASOC_ug1e2) = (xold(ASOC_ug1e2) + dt2 * P) / (1.0 + dt2 * L)


!-> ASOC_ug10

  P = YCOXY(1)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YCBENZ(1)* rct(56,k) * xnew(BENZENE) * xnew(OH)  &
    + YCTOL(1)* rct(57,k) * xnew(TOLUENE) * xnew(OH)  &
    + YCALK(1)* rct(46,k) * xnew(OH) * xnew(NC4H10)  &
    + YCOLE(1)* rct(52,k) * xnew(OH) * xnew(C3H6)  &
    + rct(82,k) * xnew(OH) * xnew(ASOC_ug1e2)
  
  L = rct(81,k) * xnew(OH)
  
  xnew(ASOC_ug10) = (xold(ASOC_ug10) + dt2 * P) / (1.0 + dt2 * L)


!-> ASOC_ug1

  P = YCOXY(0)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YCBENZ(0)* rct(56,k) * xnew(BENZENE) * xnew(OH)  &
    + YCTOL(0)* rct(57,k) * xnew(TOLUENE) * xnew(OH)  &
    + YCALK(0)* rct(46,k) * xnew(OH) * xnew(NC4H10)  &
    + YCOLE(0)* rct(52,k) * xnew(OH) * xnew(C3H6)  &
    + rct(81,k) * xnew(OH) * xnew(ASOC_ug10)
  
  L = rct(80,k) * xnew(OH)
  
  xnew(ASOC_ug1) = (xold(ASOC_ug1) + dt2 * P) / (1.0 + dt2 * L)


!-> ASOC_ng1e2

  P = rct(80,k) * xnew(OH) * xnew(ASOC_ug1)
  
  ! L = 0.0
  
  xnew(ASOC_ng1e2) = xold(ASOC_ng1e2) + dt2 * P


!-> non_C_ASOA_ug1e3

  P = YNOXY(3)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YNBENZ(3)* rct(56,k) * xnew(BENZENE) * xnew(OH)  &
    + YNTOL(3)* rct(57,k) * xnew(TOLUENE) * xnew(OH)  &
    + YNALK(3)* rct(46,k) * xnew(OH) * xnew(NC4H10)  &
    + YNOLE(3)* rct(52,k) * xnew(OH) * xnew(C3H6)
  
  L = rct(87,k) * xnew(OH)
  
  xnew(non_C_ASOA_ug1e3) = (xold(non_C_ASOA_ug1e3) + dt2 * P) / (1.0 + dt2 * L)


!-> non_C_ASOA_ug1e2

  P = YNOXY(2)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YNBENZ(2)* rct(56,k) * xnew(BENZENE) * xnew(OH)  &
    + YNTOL(2)* rct(57,k) * xnew(TOLUENE) * xnew(OH)  &
    + YNALK(2)* rct(46,k) * xnew(OH) * xnew(NC4H10)  &
    + YNOLE(2)* rct(52,k) * xnew(OH) * xnew(C3H6)  &
    + 0.9* rct(83,k) * xnew(OH) * xnew(ASOC_ug1e3)  &
    + 1.075* rct(87,k) * xnew(OH) * xnew(non_C_ASOA_ug1e3)
  
  L = rct(86,k) * xnew(OH)
  
  xnew(non_C_ASOA_ug1e2) = (xold(non_C_ASOA_ug1e2) + dt2 * P) / (1.0 + dt2 * L)


!-> non_C_ASOA_ug10

  P = YNOXY(1)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YNBENZ(1)* rct(56,k) * xnew(BENZENE) * xnew(OH)  &
    + YNTOL(1)* rct(57,k) * xnew(TOLUENE) * xnew(OH)  &
    + YNALK(1)* rct(46,k) * xnew(OH) * xnew(NC4H10)  &
    + YNOLE(1)* rct(52,k) * xnew(OH) * xnew(C3H6)  &
    + 0.9* rct(82,k) * xnew(OH) * xnew(ASOC_ug1e2)  &
    + 1.075* rct(86,k) * xnew(OH) * xnew(non_C_ASOA_ug1e2)
  
  L = rct(85,k) * xnew(OH)
  
  xnew(non_C_ASOA_ug10) = (xold(non_C_ASOA_ug10) + dt2 * P) / (1.0 + dt2 * L)


!-> non_C_ASOA_ug1

  P = YNOXY(0)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YNBENZ(0)* rct(56,k) * xnew(BENZENE) * xnew(OH)  &
    + YNTOL(0)* rct(57,k) * xnew(TOLUENE) * xnew(OH)  &
    + YNALK(0)* rct(46,k) * xnew(OH) * xnew(NC4H10)  &
    + YNOLE(0)* rct(52,k) * xnew(OH) * xnew(C3H6)  &
    + 0.9* rct(81,k) * xnew(OH) * xnew(ASOC_ug10)  &
    + 1.075* rct(85,k) * xnew(OH) * xnew(non_C_ASOA_ug10)
  
  L = rct(84,k) * xnew(OH)
  
  xnew(non_C_ASOA_ug1) = (xold(non_C_ASOA_ug1) + dt2 * P) / (1.0 + dt2 * L)


!-> non_C_ASOA_ng1e2

  P = 0.9* rct(80,k) * xnew(OH) * xnew(ASOC_ug1)  &
    + 1.075* rct(84,k) * xnew(OH) * xnew(non_C_ASOA_ug1)
  
  ! L = 0.0
  
  xnew(non_C_ASOA_ng1e2) = xold(non_C_ASOA_ng1e2) + dt2 * P


!-> SO4

  P = rcemis(SO4,k)  &
    + rct(96,k) * xnew(OH) * xnew(SO2)  &
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

  P = rct(99,k) * xnew(HNO3)  &
    + rct(100,k) * xnew(HNO3)
  
  ! L = 0.0
  
  xnew(NO3_c) = xold(NO3_c) + dt2 * P


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

if (lf_fullchem .and. k > KMAX_MID-lf_Nvert) then

!-> ASOC_ug1e3

  P_lf(1:Nd) = YCOXY(3)* 1.36e-11 * xnew_lf(1:Nd,OXYL) * xnew_lf(1:Nd,OH)  &
    + YCBENZ(3)* rct(56,k) * xnew_lf(1:Nd,BENZENE) * xnew_lf(1:Nd,OH)  &
    + YCTOL(3)* rct(57,k) * xnew_lf(1:Nd,TOLUENE) * xnew_lf(1:Nd,OH)  &
    + YCALK(3)* rct(46,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,NC4H10)  &
    + YCOLE(3)* rct(52,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,C3H6)
  
  L_lf(1:Nd) = rct(83,k) * xnew_lf(1:Nd,OH)
  
  xnew_lf(1:Nd,ASOC_ug1e3) = (xold_lf(1:Nd,ASOC_ug1e3) +  dt2 * P_lf(1:Nd)) / (1.0 + dt2 * L_lf(1:Nd))


!-> ASOC_ug1e2

  P_lf(1:Nd) = YCOXY(2)* 1.36e-11 * xnew_lf(1:Nd,OXYL) * xnew_lf(1:Nd,OH)  &
    + YCBENZ(2)* rct(56,k) * xnew_lf(1:Nd,BENZENE) * xnew_lf(1:Nd,OH)  &
    + YCTOL(2)* rct(57,k) * xnew_lf(1:Nd,TOLUENE) * xnew_lf(1:Nd,OH)  &
    + YCALK(2)* rct(46,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,NC4H10)  &
    + YCOLE(2)* rct(52,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,C3H6)  &
    + rct(83,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,ASOC_ug1e3)
  
  L_lf(1:Nd) = rct(82,k) * xnew_lf(1:Nd,OH)
  
  xnew_lf(1:Nd,ASOC_ug1e2) = (xold_lf(1:Nd,ASOC_ug1e2) +  dt2 * P_lf(1:Nd)) / (1.0 + dt2 * L_lf(1:Nd))


!-> ASOC_ug10

  P_lf(1:Nd) = YCOXY(1)* 1.36e-11 * xnew_lf(1:Nd,OXYL) * xnew_lf(1:Nd,OH)  &
    + YCBENZ(1)* rct(56,k) * xnew_lf(1:Nd,BENZENE) * xnew_lf(1:Nd,OH)  &
    + YCTOL(1)* rct(57,k) * xnew_lf(1:Nd,TOLUENE) * xnew_lf(1:Nd,OH)  &
    + YCALK(1)* rct(46,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,NC4H10)  &
    + YCOLE(1)* rct(52,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,C3H6)  &
    + rct(82,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,ASOC_ug1e2)
  
  L_lf(1:Nd) = rct(81,k) * xnew_lf(1:Nd,OH)
  
  xnew_lf(1:Nd,ASOC_ug10) = (xold_lf(1:Nd,ASOC_ug10) +  dt2 * P_lf(1:Nd)) / (1.0 + dt2 * L_lf(1:Nd))


!-> ASOC_ug1

  P_lf(1:Nd) = YCOXY(0)* 1.36e-11 * xnew_lf(1:Nd,OXYL) * xnew_lf(1:Nd,OH)  &
    + YCBENZ(0)* rct(56,k) * xnew_lf(1:Nd,BENZENE) * xnew_lf(1:Nd,OH)  &
    + YCTOL(0)* rct(57,k) * xnew_lf(1:Nd,TOLUENE) * xnew_lf(1:Nd,OH)  &
    + YCALK(0)* rct(46,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,NC4H10)  &
    + YCOLE(0)* rct(52,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,C3H6)  &
    + rct(81,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,ASOC_ug10)
  
  L_lf(1:Nd) = rct(80,k) * xnew_lf(1:Nd,OH)
  
  xnew_lf(1:Nd,ASOC_ug1) = (xold_lf(1:Nd,ASOC_ug1) +  dt2 * P_lf(1:Nd)) / (1.0 + dt2 * L_lf(1:Nd))


!-> ASOC_ng1e2

  P_lf(1:Nd) = rct(80,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,ASOC_ug1)
  
  ! L_lf(1:Nd) = 0.0
  
  xnew_lf(1:Nd,ASOC_ng1e2) = xold_lf(1:Nd,ASOC_ng1e2) +  dt2 * P_lf(1:Nd)


!-> non_C_ASOA_ug1e3

  P_lf(1:Nd) = YNOXY(3)* 1.36e-11 * xnew_lf(1:Nd,OXYL) * xnew_lf(1:Nd,OH)  &
    + YNBENZ(3)* rct(56,k) * xnew_lf(1:Nd,BENZENE) * xnew_lf(1:Nd,OH)  &
    + YNTOL(3)* rct(57,k) * xnew_lf(1:Nd,TOLUENE) * xnew_lf(1:Nd,OH)  &
    + YNALK(3)* rct(46,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,NC4H10)  &
    + YNOLE(3)* rct(52,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,C3H6)
  
  L_lf(1:Nd) = rct(87,k) * xnew_lf(1:Nd,OH)
  
  xnew_lf(1:Nd,non_C_ASOA_ug1e3) = (xold_lf(1:Nd,non_C_ASOA_ug1e3) +  dt2 * P_lf(1:Nd)) / (1.0 + dt2 * L_lf(1:Nd))


!-> non_C_ASOA_ug1e2

  P_lf(1:Nd) = YNOXY(2)* 1.36e-11 * xnew_lf(1:Nd,OXYL) * xnew_lf(1:Nd,OH)  &
    + YNBENZ(2)* rct(56,k) * xnew_lf(1:Nd,BENZENE) * xnew_lf(1:Nd,OH)  &
    + YNTOL(2)* rct(57,k) * xnew_lf(1:Nd,TOLUENE) * xnew_lf(1:Nd,OH)  &
    + YNALK(2)* rct(46,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,NC4H10)  &
    + YNOLE(2)* rct(52,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,C3H6)  &
    + 0.9* rct(83,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,ASOC_ug1e3)  &
    + 1.075* rct(87,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,non_C_ASOA_ug1e3)
  
  L_lf(1:Nd) = rct(86,k) * xnew_lf(1:Nd,OH)
  
  xnew_lf(1:Nd,non_C_ASOA_ug1e2) = (xold_lf(1:Nd,non_C_ASOA_ug1e2) +  dt2 * P_lf(1:Nd)) / (1.0 + dt2 * L_lf(1:Nd))


!-> non_C_ASOA_ug10

  P_lf(1:Nd) = YNOXY(1)* 1.36e-11 * xnew_lf(1:Nd,OXYL) * xnew_lf(1:Nd,OH)  &
    + YNBENZ(1)* rct(56,k) * xnew_lf(1:Nd,BENZENE) * xnew_lf(1:Nd,OH)  &
    + YNTOL(1)* rct(57,k) * xnew_lf(1:Nd,TOLUENE) * xnew_lf(1:Nd,OH)  &
    + YNALK(1)* rct(46,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,NC4H10)  &
    + YNOLE(1)* rct(52,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,C3H6)  &
    + 0.9* rct(82,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,ASOC_ug1e2)  &
    + 1.075* rct(86,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,non_C_ASOA_ug1e2)
  
  L_lf(1:Nd) = rct(85,k) * xnew_lf(1:Nd,OH)
  
  xnew_lf(1:Nd,non_C_ASOA_ug10) = (xold_lf(1:Nd,non_C_ASOA_ug10) +  dt2 * P_lf(1:Nd)) / (1.0 + dt2 * L_lf(1:Nd))


!-> non_C_ASOA_ug1

  P_lf(1:Nd) = YNOXY(0)* 1.36e-11 * xnew_lf(1:Nd,OXYL) * xnew_lf(1:Nd,OH)  &
    + YNBENZ(0)* rct(56,k) * xnew_lf(1:Nd,BENZENE) * xnew_lf(1:Nd,OH)  &
    + YNTOL(0)* rct(57,k) * xnew_lf(1:Nd,TOLUENE) * xnew_lf(1:Nd,OH)  &
    + YNALK(0)* rct(46,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,NC4H10)  &
    + YNOLE(0)* rct(52,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,C3H6)  &
    + 0.9* rct(81,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,ASOC_ug10)  &
    + 1.075* rct(85,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,non_C_ASOA_ug10)
  
  L_lf(1:Nd) = rct(84,k) * xnew_lf(1:Nd,OH)
  
  xnew_lf(1:Nd,non_C_ASOA_ug1) = (xold_lf(1:Nd,non_C_ASOA_ug1) +  dt2 * P_lf(1:Nd)) / (1.0 + dt2 * L_lf(1:Nd))


!-> non_C_ASOA_ng1e2

  P_lf(1:Nd) = 0.9* rct(80,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,ASOC_ug1)  &
    + 1.075* rct(84,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,non_C_ASOA_ug1)
  
  ! L_lf(1:Nd) = 0.0
  
  xnew_lf(1:Nd,non_C_ASOA_ng1e2) = xold_lf(1:Nd,non_C_ASOA_ng1e2) +  dt2 * P_lf(1:Nd)


!-> SO4

  P_lf(1:Nd) = rcemis(SO4,k)  &
    + rct(96,k)*fgasso2_lf(1:Nd,k) * xnew_lf(1:Nd,OH) * xnew_lf(1:Nd,SO2)  &
    + AQRCK_LF(1:ND,ICLRC1,K) * xnew_lf(1:Nd,SO2) * xnew_lf(1:Nd,H2O2)  &
    + AQRCK_LF(1:ND,ICLRC2,K) * xnew_lf(1:Nd,SO2) * xnew_lf(1:Nd,O3)  &
    + AQRCK_LF(1:ND,ICLRC3,K) * xnew_lf(1:Nd,SO2)
  
  ! L_lf(1:Nd) = 0.0
  
P_lf(NSPEC_deriv_lf+1:Nd) = P_lf(NSPEC_deriv_lf+1:Nd) + rcemis_lf(1:N_lf_derivemis,SO4) 

  xnew_lf(1:Nd,SO4) = xold_lf(1:Nd,SO4) +  dt2 * P_lf(1:Nd)


!-> NH3

  P_lf(1:Nd) = rcemis(NH3,k)
  
  ! L_lf(1:Nd) = 0.0
  
P_lf(NSPEC_deriv_lf+1:Nd) = P_lf(NSPEC_deriv_lf+1:Nd) + rcemis_lf(1:N_lf_derivemis,NH3) 

  xnew_lf(1:Nd,NH3) = xold_lf(1:Nd,NH3) +  dt2 * P_lf(1:Nd)


!-> NO3_f

  ! P_lf(1:Nd) = 0.0
  
  ! L_lf(1:Nd) = 0.0
  
  ! Nothing to do for NO3_f! xnew_lf(1:Nd,NO3_f) = max(0.0, xold_lf(1:Nd,NO3_f))


!-> NO3_c

  P_lf(1:Nd) = rct(99,k) * xnew_lf(1:Nd,HNO3)  &
    + rct(100,k) * xnew_lf(1:Nd,HNO3)
  
  ! L_lf(1:Nd) = 0.0
  
  xnew_lf(1:Nd,NO3_c) = xold_lf(1:Nd,NO3_c) +  dt2 * P_lf(1:Nd)

!-> ASOC_ug1e3

  P_lf(Nd+1:Nd+NSOA) = YCOXY(3)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YCBENZ(3)* rct(56,k) * xnew(BENZENE) * xnew(OH)  &
    + YCTOL(3)* rct(57,k) * xnew(TOLUENE) * xnew(OH)  &
    + YCALK(3)* rct(46,k) * xnew(OH) * xnew(NC4H10)  &
    + YCOLE(3)* rct(52,k) * xnew(OH) * xnew(C3H6)
  
  L_lf(Nd+1:Nd+NSOA) = rct(83,k) * xnew(OH)
  
  xnew_lf(Nd+1:Nd+NSOA,ASOC_ug1e3) = (xold_lf(Nd+1:Nd+NSOA,ASOC_ug1e3) + dt2 * P_lf(Nd+1:Nd+NSOA)) / (1.0 + dt2 * L_lf(Nd+1:Nd+NSOA))


!-> ASOC_ug1e2

  P_lf(Nd+1:Nd+NSOA) = YCOXY(2)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YCBENZ(2)* rct(56,k) * xnew(BENZENE) * xnew(OH)  &
    + YCTOL(2)* rct(57,k) * xnew(TOLUENE) * xnew(OH)  &
    + YCALK(2)* rct(46,k) * xnew(OH) * xnew(NC4H10)  &
    + YCOLE(2)* rct(52,k) * xnew(OH) * xnew(C3H6)  &
    + rct(83,k) * xnew(OH) * xnew_lf(Nd+1:Nd+NSOA,ASOC_ug1e3)
  
  L_lf(Nd+1:Nd+NSOA) = rct(82,k) * xnew(OH)
  
  xnew_lf(Nd+1:Nd+NSOA,ASOC_ug1e2) = (xold_lf(Nd+1:Nd+NSOA,ASOC_ug1e2) + dt2 * P_lf(Nd+1:Nd+NSOA)) / (1.0 + dt2 * L_lf(Nd+1:Nd+NSOA))


!-> ASOC_ug10

  P_lf(Nd+1:Nd+NSOA) = YCOXY(1)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YCBENZ(1)* rct(56,k) * xnew(BENZENE) * xnew(OH)  &
    + YCTOL(1)* rct(57,k) * xnew(TOLUENE) * xnew(OH)  &
    + YCALK(1)* rct(46,k) * xnew(OH) * xnew(NC4H10)  &
    + YCOLE(1)* rct(52,k) * xnew(OH) * xnew(C3H6)  &
    + rct(82,k) * xnew(OH) * xnew_lf(Nd+1:Nd+NSOA,ASOC_ug1e2)
  
  L_lf(Nd+1:Nd+NSOA) = rct(81,k) * xnew(OH)
  
  xnew_lf(Nd+1:Nd+NSOA,ASOC_ug10) = (xold_lf(Nd+1:Nd+NSOA,ASOC_ug10) + dt2 * P_lf(Nd+1:Nd+NSOA)) / (1.0 + dt2 * L_lf(Nd+1:Nd+NSOA))


!-> ASOC_ug1

  P_lf(Nd+1:Nd+NSOA) = YCOXY(0)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YCBENZ(0)* rct(56,k) * xnew(BENZENE) * xnew(OH)  &
    + YCTOL(0)* rct(57,k) * xnew(TOLUENE) * xnew(OH)  &
    + YCALK(0)* rct(46,k) * xnew(OH) * xnew(NC4H10)  &
    + YCOLE(0)* rct(52,k) * xnew(OH) * xnew(C3H6)  &
    + rct(81,k) * xnew(OH) * xnew_lf(Nd+1:Nd+NSOA,ASOC_ug10)
  
  L_lf(Nd+1:Nd+NSOA) = rct(80,k) * xnew(OH)
  
  xnew_lf(Nd+1:Nd+NSOA,ASOC_ug1) = (xold_lf(Nd+1:Nd+NSOA,ASOC_ug1) + dt2 * P_lf(Nd+1:Nd+NSOA)) / (1.0 + dt2 * L_lf(Nd+1:Nd+NSOA))


!-> ASOC_ng1e2

  P_lf(Nd+1:Nd+NSOA) = rct(80,k) * xnew(OH) * xnew_lf(Nd+1:Nd+NSOA,ASOC_ug1)
  
  ! L_lf(Nd+1:Nd+NSOA) = 0.0
  
  xnew_lf(Nd+1:Nd+NSOA,ASOC_ng1e2) = xold_lf(Nd+1:Nd+NSOA,ASOC_ng1e2) + dt2 * P_lf(Nd+1:Nd+NSOA)


!-> non_C_ASOA_ug1e3

  P_lf(Nd+1:Nd+NSOA) = YNOXY(3)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YNBENZ(3)* rct(56,k) * xnew(BENZENE) * xnew(OH)  &
    + YNTOL(3)* rct(57,k) * xnew(TOLUENE) * xnew(OH)  &
    + YNALK(3)* rct(46,k) * xnew(OH) * xnew(NC4H10)  &
    + YNOLE(3)* rct(52,k) * xnew(OH) * xnew(C3H6)
  
  L_lf(Nd+1:Nd+NSOA) = rct(87,k) * xnew(OH)
  
  xnew_lf(Nd+1:Nd+NSOA,non_C_ASOA_ug1e3) = (xold_lf(Nd+1:Nd+NSOA,non_C_ASOA_ug1e3) + dt2 * P_lf(Nd+1:Nd+NSOA)) / (1.0 + dt2 * L_lf(Nd+1:Nd+NSOA))


!-> non_C_ASOA_ug1e2

  P_lf(Nd+1:Nd+NSOA) = YNOXY(2)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YNBENZ(2)* rct(56,k) * xnew(BENZENE) * xnew(OH)  &
    + YNTOL(2)* rct(57,k) * xnew(TOLUENE) * xnew(OH)  &
    + YNALK(2)* rct(46,k) * xnew(OH) * xnew(NC4H10)  &
    + YNOLE(2)* rct(52,k) * xnew(OH) * xnew(C3H6)  &
    + 0.9* rct(83,k) * xnew(OH) * xnew_lf(Nd+1:Nd+NSOA,ASOC_ug1e3)  &
    + 1.075* rct(87,k) * xnew(OH) * xnew_lf(Nd+1:Nd+NSOA,non_C_ASOA_ug1e3)
  
  L_lf(Nd+1:Nd+NSOA) = rct(86,k) * xnew(OH)
  
  xnew_lf(Nd+1:Nd+NSOA,non_C_ASOA_ug1e2) = (xold_lf(Nd+1:Nd+NSOA,non_C_ASOA_ug1e2) + dt2 * P_lf(Nd+1:Nd+NSOA)) / (1.0 + dt2 * L_lf(Nd+1:Nd+NSOA))


!-> non_C_ASOA_ug10

  P_lf(Nd+1:Nd+NSOA) = YNOXY(1)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YNBENZ(1)* rct(56,k) * xnew(BENZENE) * xnew(OH)  &
    + YNTOL(1)* rct(57,k) * xnew(TOLUENE) * xnew(OH)  &
    + YNALK(1)* rct(46,k) * xnew(OH) * xnew(NC4H10)  &
    + YNOLE(1)* rct(52,k) * xnew(OH) * xnew(C3H6)  &
    + 0.9* rct(82,k) * xnew(OH) * xnew_lf(Nd+1:Nd+NSOA,ASOC_ug1e2)  &
    + 1.075* rct(86,k) * xnew(OH) * xnew_lf(Nd+1:Nd+NSOA,non_C_ASOA_ug1e2)
  
  L_lf(Nd+1:Nd+NSOA) = rct(85,k) * xnew(OH)
  
  xnew_lf(Nd+1:Nd+NSOA,non_C_ASOA_ug10) = (xold_lf(Nd+1:Nd+NSOA,non_C_ASOA_ug10) + dt2 * P_lf(Nd+1:Nd+NSOA)) / (1.0 + dt2 * L_lf(Nd+1:Nd+NSOA))


!-> non_C_ASOA_ug1

  P_lf(Nd+1:Nd+NSOA) = YNOXY(0)* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + YNBENZ(0)* rct(56,k) * xnew(BENZENE) * xnew(OH)  &
    + YNTOL(0)* rct(57,k) * xnew(TOLUENE) * xnew(OH)  &
    + YNALK(0)* rct(46,k) * xnew(OH) * xnew(NC4H10)  &
    + YNOLE(0)* rct(52,k) * xnew(OH) * xnew(C3H6)  &
    + 0.9* rct(81,k) * xnew(OH) * xnew_lf(Nd+1:Nd+NSOA,ASOC_ug10)  &
    + 1.075* rct(85,k) * xnew(OH) * xnew_lf(Nd+1:Nd+NSOA,non_C_ASOA_ug10)
  
  L_lf(Nd+1:Nd+NSOA) = rct(84,k) * xnew(OH)
  
  xnew_lf(Nd+1:Nd+NSOA,non_C_ASOA_ug1) = (xold_lf(Nd+1:Nd+NSOA,non_C_ASOA_ug1) + dt2 * P_lf(Nd+1:Nd+NSOA)) / (1.0 + dt2 * L_lf(Nd+1:Nd+NSOA))


!-> non_C_ASOA_ng1e2

  P_lf(Nd+1:Nd+NSOA) = 0.9* rct(80,k) * xnew(OH) * xnew_lf(Nd+1:Nd+NSOA,ASOC_ug1)  &
    + 1.075* rct(84,k) * xnew(OH) * xnew_lf(Nd+1:Nd+NSOA,non_C_ASOA_ug1)
  
  ! L_lf(Nd+1:Nd+NSOA) = 0.0
  
  xnew_lf(Nd+1:Nd+NSOA,non_C_ASOA_ng1e2) = xold_lf(Nd+1:Nd+NSOA,non_C_ASOA_ng1e2) + dt2 * P_lf(Nd+1:Nd+NSOA)


endif
!-> BSOC_ug1e3

  P = YCISOP(3)* rct(64,k) * xnew(OH) * xnew(C5H8)  &
    + YCTERP(3)* rct(75,k) * xnew(OH) * xnew(APINENE)  &
    + YCTERP(3)* rct(78,k) * xnew(O3) * xnew(APINENE)  &
    + YCTERP(3)* rct(79,k) * xnew(NO3) * xnew(APINENE)
  
  L = rct(91,k) * xnew(OH)
  
  xnew(BSOC_ug1e3) = (xold(BSOC_ug1e3) + dt2 * P) / (1.0 + dt2 * L)


!-> BSOC_ug1e2

  P = YCISOP(2)* rct(64,k) * xnew(OH) * xnew(C5H8)  &
    + YCTERP(2)* rct(75,k) * xnew(OH) * xnew(APINENE)  &
    + YCTERP(2)* rct(78,k) * xnew(O3) * xnew(APINENE)  &
    + YCTERP(2)* rct(79,k) * xnew(NO3) * xnew(APINENE)  &
    + rct(91,k) * xnew(OH) * xnew(BSOC_ug1e3)
  
  L = rct(90,k) * xnew(OH)
  
  xnew(BSOC_ug1e2) = (xold(BSOC_ug1e2) + dt2 * P) / (1.0 + dt2 * L)


!-> BSOC_ug10

  P = YCISOP(1)* rct(64,k) * xnew(OH) * xnew(C5H8)  &
    + YCTERP(1)* rct(75,k) * xnew(OH) * xnew(APINENE)  &
    + YCTERP(1)* rct(78,k) * xnew(O3) * xnew(APINENE)  &
    + YCTERP(1)* rct(79,k) * xnew(NO3) * xnew(APINENE)  &
    + rct(90,k) * xnew(OH) * xnew(BSOC_ug1e2)
  
  L = rct(89,k) * xnew(OH)
  
  xnew(BSOC_ug10) = (xold(BSOC_ug10) + dt2 * P) / (1.0 + dt2 * L)


!-> BSOC_ug1

  P = YCISOP(0)* rct(64,k) * xnew(OH) * xnew(C5H8)  &
    + YCTERP(0)* rct(75,k) * xnew(OH) * xnew(APINENE)  &
    + YCTERP(0)* rct(78,k) * xnew(O3) * xnew(APINENE)  &
    + YCTERP(0)* rct(79,k) * xnew(NO3) * xnew(APINENE)  &
    + rct(89,k) * xnew(OH) * xnew(BSOC_ug10)
  
  L = rct(88,k) * xnew(OH)
  
  xnew(BSOC_ug1) = (xold(BSOC_ug1) + dt2 * P) / (1.0 + dt2 * L)


!-> BSOC_ng1e2

  P = rct(88,k) * xnew(OH) * xnew(BSOC_ug1)
  
  ! L = 0.0
  
  xnew(BSOC_ng1e2) = xold(BSOC_ng1e2) + dt2 * P


!-> non_C_BSOA_ug1e3

  P = YNISOP(3)* rct(64,k) * xnew(OH) * xnew(C5H8)  &
    + YNTERP(3)* rct(75,k) * xnew(OH) * xnew(APINENE)  &
    + YNTERP(3)* rct(78,k) * xnew(O3) * xnew(APINENE)  &
    + YNTERP(3)* rct(79,k) * xnew(NO3) * xnew(APINENE)
  
  L = rct(95,k) * xnew(OH)
  
  xnew(non_C_BSOA_ug1e3) = (xold(non_C_BSOA_ug1e3) + dt2 * P) / (1.0 + dt2 * L)


!-> non_C_BSOA_ug1e2

  P = YNISOP(2)* rct(64,k) * xnew(OH) * xnew(C5H8)  &
    + YNTERP(2)* rct(75,k) * xnew(OH) * xnew(APINENE)  &
    + YNTERP(2)* rct(78,k) * xnew(O3) * xnew(APINENE)  &
    + YNTERP(2)* rct(79,k) * xnew(NO3) * xnew(APINENE)  &
    + 0.9* rct(91,k) * xnew(OH) * xnew(BSOC_ug1e3)  &
    + 1.075* rct(95,k) * xnew(OH) * xnew(non_C_BSOA_ug1e3)
  
  L = rct(94,k) * xnew(OH)
  
  xnew(non_C_BSOA_ug1e2) = (xold(non_C_BSOA_ug1e2) + dt2 * P) / (1.0 + dt2 * L)


!-> non_C_BSOA_ug10

  P = YNISOP(1)* rct(64,k) * xnew(OH) * xnew(C5H8)  &
    + YNTERP(1)* rct(75,k) * xnew(OH) * xnew(APINENE)  &
    + YNTERP(1)* rct(78,k) * xnew(O3) * xnew(APINENE)  &
    + YNTERP(1)* rct(79,k) * xnew(NO3) * xnew(APINENE)  &
    + 0.9* rct(90,k) * xnew(OH) * xnew(BSOC_ug1e2)  &
    + 1.075* rct(94,k) * xnew(OH) * xnew(non_C_BSOA_ug1e2)
  
  L = rct(93,k) * xnew(OH)
  
  xnew(non_C_BSOA_ug10) = (xold(non_C_BSOA_ug10) + dt2 * P) / (1.0 + dt2 * L)


!-> non_C_BSOA_ug1

  P = YNISOP(0)* rct(64,k) * xnew(OH) * xnew(C5H8)  &
    + YNTERP(0)* rct(75,k) * xnew(OH) * xnew(APINENE)  &
    + YNTERP(0)* rct(78,k) * xnew(O3) * xnew(APINENE)  &
    + YNTERP(0)* rct(79,k) * xnew(NO3) * xnew(APINENE)  &
    + 0.9* rct(89,k) * xnew(OH) * xnew(BSOC_ug10)  &
    + 1.075* rct(93,k) * xnew(OH) * xnew(non_C_BSOA_ug10)
  
  L = rct(92,k) * xnew(OH)
  
  xnew(non_C_BSOA_ug1) = (xold(non_C_BSOA_ug1) + dt2 * P) / (1.0 + dt2 * L)


!-> non_C_BSOA_ng1e2

  P = 0.9* rct(88,k) * xnew(OH) * xnew(BSOC_ug1)  &
    + 1.075* rct(92,k) * xnew(OH) * xnew(non_C_BSOA_ug1)
  
  ! L = 0.0
  
  xnew(non_C_BSOA_ng1e2) = xold(non_C_BSOA_ng1e2) + dt2 * P


!-> NH4_f

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for NH4_f! xnew(NH4_f) = max(0.0, xold(NH4_f))


!-> OM25_bgnd

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for OM25_bgnd! xnew(OM25_bgnd) = max(0.0, xold(OM25_bgnd))


!-> OM25_p

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for OM25_p! xnew(OM25_p) = max(0.0, xold(OM25_p))


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


!-> ffire_c

  P = rcemis(FFIRE_C,k)
  
  ! L = 0.0
  
  xnew(ffire_c) = xold(ffire_c) + dt2 * P


!-> SeaSalt_f

  P = rcemis(SEASALT_F,K)
  
  ! L = 0.0
  
  xnew(SeaSalt_f) = xold(SeaSalt_f) + dt2 * P


!-> SeaSalt_c

  P = rcemis(SEASALT_C,K)
  
  ! L = 0.0
  
  xnew(SeaSalt_c) = xold(SeaSalt_c) + dt2 * P


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


!-> Ash_f

  P = rcemis(ASH_F,K)
  
  ! L = 0.0
  
  xnew(Ash_f) = xold(Ash_f) + dt2 * P


!-> Ash_c

  P = rcemis(ASH_C,K)
  
  ! L = 0.0
  
  xnew(Ash_c) = xold(Ash_c) + dt2 * P


!-> POM_f_Res

  P = rcemis(POM_F_RES,k)
  
  ! L = 0.0
  
  xnew(POM_f_Res) = xold(POM_f_Res) + dt2 * P


!-> POM_c_Res

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for POM_c_Res! xnew(POM_c_Res) = max(0.0, xold(POM_c_Res))


!-> POM_f_nonRes

  P = rcemis(POM_F_NONRES,k)
  
  ! L = 0.0
  
  xnew(POM_f_nonRes) = xold(POM_f_nonRes) + dt2 * P


!-> POM_c_nonRes

  P = rcemis(POM_C_NONRES,k)
  
  ! L = 0.0
  
  xnew(POM_c_nonRes) = xold(POM_c_nonRes) + dt2 * P


!-> EC_f_Res_new

  P = rcemis(EC_F_RES_NEW,k)
  
  L = rct(102,k)
  
  xnew(EC_f_Res_new) = (xold(EC_f_Res_new) + dt2 * P) / (1.0 + dt2 * L)


!-> EC_f_Res_age

  P = rcemis(EC_F_RES_AGE,k)  &
    + rct(102,k) * xnew(EC_f_Res_new)
  
  ! L = 0.0
  
  xnew(EC_f_Res_age) = xold(EC_f_Res_age) + dt2 * P


!-> EC_c_Res

  P = rcemis(EC_C_RES,k)
  
  ! L = 0.0
  
  xnew(EC_c_Res) = xold(EC_c_Res) + dt2 * P


!-> EC_f_nonRes_new

  P = rcemis(EC_F_NONRES_NEW,k)
  
  L = rct(102,k)
  
  xnew(EC_f_nonRes_new) = (xold(EC_f_nonRes_new) + dt2 * P) / (1.0 + dt2 * L)


!-> EC_f_nonRes_age

  P = rcemis(EC_F_NONRES_AGE,k)  &
    + rct(102,k) * xnew(EC_f_nonRes_new)
  
  ! L = 0.0
  
  xnew(EC_f_nonRes_age) = xold(EC_f_nonRes_age) + dt2 * P


!-> EC_c_nonRes

  P = rcemis(EC_C_NONRES,k)
  
  ! L = 0.0
  
  xnew(EC_c_nonRes) = xold(EC_c_nonRes) + dt2 * P


!-> remPPM25_nonRes

  P = rcemis(REMPPM25_NONRES,k)
  
  ! L = 0.0
  
  xnew(remPPM25_nonRes) = xold(remPPM25_nonRes) + dt2 * P


!-> remPPM25_Res

  P = rcemis(REMPPM25_RES,k)
  
  ! L = 0.0
  
  xnew(remPPM25_Res) = xold(remPPM25_Res) + dt2 * P


!-> remPPM_c_nonRes

  P = rcemis(REMPPM_C_NONRES,k)
  
  ! L = 0.0
  
  xnew(remPPM_c_nonRes) = xold(remPPM_c_nonRes) + dt2 * P


!-> remPPM_c_Res

  P = rcemis(REMPPM_C_RES,k)
  
  ! L = 0.0
  
  xnew(remPPM_c_Res) = xold(remPPM_c_Res) + dt2 * P


!-> SQT_SOA_NV

  P = 0.00383* RCBIO(2,K)
  
  ! L = 0.0
  
  xnew(SQT_SOA_NV) = xold(SQT_SOA_NV) + dt2 * P


