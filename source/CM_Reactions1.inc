
!-> OD 
      ! P = 0.0
 

      L =  &
          rct(2,k) 

      xnew(OD)=  xold(OD) / ( 1.0 + dt2 * L )

!-> OP 
      ! P = 0.0
 

      L =  &
          rct(1,k) 

      xnew(OP)=  xold(OP) / ( 1.0 + dt2 * L )

!-> OH 

      P =   &
          2.*rct(2,k) * xnew(OD )  &
        + rct(6,k) * xnew(NO ) * xnew(HO2 ) 

      L =  &
          rct(4,k)* xnew(VOC )  &
        + rct(5,k)* xnew(CH3CHO )  &
        + 1.11e-11* xnew(NO2 ) 

      xnew(OH)=  ( xold(OH) + dt2 * P) /(1.0 + dt2*L )

!-> HO2 

      P =   &
          rct(7,k) * xnew(NO ) * xnew(RO2 ) 

      L =  &
          rct(6,k)* xnew(NO ) 

      xnew(HO2)=  ( xold(HO2) + dt2 * P) /(1.0 + dt2*L )

!-> RO2 

      P =   &
          rct(4,k) * xnew(OH ) * xnew(VOC )  &
        + rct(8,k) * xnew(CH3COO2 ) * xnew(NO ) 

      L =  &
          rct(7,k)* xnew(NO ) 

      xnew(RO2)=  ( xold(RO2) + dt2 * P) /(1.0 + dt2*L )

!-> TRACER1 

      P =   &
          rcemis(TRACER1,k) 
      ! L = 0.0
 

      xnew(TRACER1)=  xold(TRACER1) + dt2 * P

!-> TRACER2 

      P =   &
          rcemis(TRACER2,k) 
      ! L = 0.0
 

      xnew(TRACER2)=  xold(TRACER2) + dt2 * P

!-> O3 

      P =   &
          rct(1,k) * xnew(OP ) 

      L =  &
          rct(3,k)* xnew(NO ) 

      xnew(O3)=  ( xold(O3) + dt2 * P) /(1.0 + dt2*L )

!-> NO 

      P =   &
          rcemis(NO,k)  &
        + 0 !Skip bio rate since rcemis exists 

      L =  &
          rct(3,k)* xnew(O3 )  &
        + rct(6,k)* xnew(HO2 )  &
        + rct(7,k)* xnew(RO2 )  &
        + rct(8,k)* xnew(CH3COO2 ) 

      xnew(NO)=  ( xold(NO) + dt2 * P) /(1.0 + dt2*L )

!-> NO2 

      P =   &
          rct(3,k) * xnew(O3 ) * xnew(NO )  &
        + rct(6,k) * xnew(NO ) * xnew(HO2 )  &
        + rct(7,k) * xnew(NO ) * xnew(RO2 )  &
        + rct(8,k) * xnew(CH3COO2 ) * xnew(NO )  &
        + rct(9,k) * xnew(PAN )  &
        + rcemis(NO2,k) 

      L =  &
          1.11e-11* xnew(OH )  &
        + 4.7e-12* xnew(CH3COO2 ) 

      xnew(NO2)=  ( xold(NO2) + dt2 * P) /(1.0 + dt2*L )

!-> SHIPNOX 

      P =   &
          rcemis(SHIPNOX,k) 
      ! L = 0.0
 

      xnew(SHIPNOX)=  xold(SHIPNOX) + dt2 * P

!-> HNO3 

      P =   &
          1.11e-11 * xnew(NO2 ) * xnew(OH ) 
      ! L = 0.0
 

      xnew(HNO3)=  xold(HNO3) + dt2 * P

!-> VOC 
      ! P = 0.0
 

      L =  &
          rct(4,k)* xnew(OH ) 

      xnew(VOC)=  xold(VOC) / ( 1.0 + dt2 * L )

!-> CH3CHO 

      P =   &
          rct(7,k) * xnew(NO ) * xnew(RO2 )  &
        + rcemis(CH3CHO,k) 

      L =  &
          rct(5,k)* xnew(OH ) 

      xnew(CH3CHO)=  ( xold(CH3CHO) + dt2 * P) /(1.0 + dt2*L )

!-> PAN 

      P =   &
          4.7e-12 * xnew(CH3COO2 ) * xnew(NO2 ) 

      L =  &
          rct(9,k) 

      xnew(PAN)=  ( xold(PAN) + dt2 * P) /(1.0 + dt2*L )

!-> CH3COO2 

      P =   &
          rct(5,k) * xnew(OH ) * xnew(CH3CHO )  &
        + rct(9,k) * xnew(PAN ) 

      L =  &
          rct(8,k)* xnew(NO )  &
        + 4.7e-12* xnew(NO2 ) 

      xnew(CH3COO2)=  ( xold(CH3COO2) + dt2 * P) /(1.0 + dt2*L )

!-> SO2 

      P =   &
          rcemis(SO2,k) 
      ! L = 0.0
 

      xnew(SO2)=  xold(SO2) + dt2 * P

!-> SO4 

      P =   &
          rcemis(SO4,k) 
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
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for NO3_C! xnew(NO3_C)=  max(0.0, xold(NO3_C))

!-> NH4_F 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for NH4_F! xnew(NH4_F)=  max(0.0, xold(NH4_F))
