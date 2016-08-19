
!-> SO4 

      P =   &
          2e-12*AQRCK(ICLOHSO2,K) * xnew(OH ) * xnew(SO2 )  &
        + AQRCK(ICLRC1,K) * xnew(SO2 ) * xnew(H2O2 )  &
        + AQRCK(ICLRC2,K) * xnew(SO2 ) * xnew(O3 )  &
        + AQRCK(ICLRC3,K) * xnew(SO2 )  &
        + rcemis(SO4,k) 
      ! L = 0.0
 

      xnew(SO4)=  max(0.0,  xold(SO4) + dt2 * P)

!-> NH3 

      P =   &
          rcemis(NH3,k) 
      ! L = 0.0
 

      xnew(NH3)=  max(0.0,  xold(NH3) + dt2 * P)

!-> NO3_F 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for NO3_F! xnew(NO3_F)=  max(0.0, xold(NO3_F))

!-> NO3_C 

      P =   &
          rcmisc(19,k) * xnew(HNO3 ) 
      ! L = 0.0
 

      xnew(NO3_C)=  max(0.0,  xold(NO3_C) + dt2 * P)

!-> NH4_F 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for NH4_F! xnew(NH4_F)=  max(0.0, xold(NH4_F))

!-> PPM25 

      P =   &
          rcemis(PPM25,k) 
      ! L = 0.0
 

      xnew(PPM25)=  max(0.0,  xold(PPM25) + dt2 * P)

!-> PPM25_FIRE 

      P =   &
          rcemis(PPM25_FIRE,k) 
      ! L = 0.0
 

      xnew(PPM25_FIRE)=  max(0.0,  xold(PPM25_FIRE) + dt2 * P)

!-> PPM_C 

      P =   &
          rcemis(PPM_C,k) 
      ! L = 0.0
 

      xnew(PPM_C)=  max(0.0,  xold(PPM_C) + dt2 * P)

!-> SEASALT_F 

      P =   &
          RCSS(QSSFI,K) 
      ! L = 0.0
 

      xnew(SEASALT_F)=  max(0.0,  xold(SEASALT_F) + dt2 * P)

!-> SEASALT_C 

      P =   &
          RCSS(QSSCO,K) 
      ! L = 0.0
 

      xnew(SEASALT_C)=  max(0.0,  xold(SEASALT_C) + dt2 * P)

!-> SEASALT_G 

      P =   &
          RCSS(QSSGI,K) 
      ! L = 0.0
 

      xnew(SEASALT_G)=  max(0.0,  xold(SEASALT_G) + dt2 * P)

!-> DUST_NAT_F 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for DUST_NAT_F! xnew(DUST_NAT_F)=  max(0.0, xold(DUST_NAT_F))

!-> DUST_NAT_C 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for DUST_NAT_C! xnew(DUST_NAT_C)=  max(0.0, xold(DUST_NAT_C))

!-> RN222 

      P =   &
          RC_RN222(K) 

      L =  &
          2.1e-6 

      xnew(RN222)=  max(0.0, ( xold(RN222) + dt2 * P)) /(1.0 + dt2*L )

!-> PB210 

      P =   &
          2.1e-6 * xnew(RN222 ) 

      L =  &
          rcmisc(18,k) 

      xnew(PB210)=  max(0.0, ( xold(PB210) + dt2 * P)) /(1.0 + dt2*L )
