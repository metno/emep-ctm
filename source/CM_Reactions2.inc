
!-> CO 

      P =   &
          rcemis(CO,k) 
      ! L = 0.0
 

      xnew(CO)=  xold(CO) + dt2 * P

!-> HCHO 

      P =   &
          rcemis(HCHO,k) 
      ! L = 0.0
 

      xnew(HCHO)=  xold(HCHO) + dt2 * P

!-> CH4 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for CH4! xnew(CH4)=  max(0.0, xold(CH4))

!-> H2 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for H2! xnew(H2)=  max(0.0, xold(H2))

!-> C2H6 

      P =   &
          rcemis(C2H6,k) 
      ! L = 0.0
 

      xnew(C2H6)=  xold(C2H6) + dt2 * P

!-> NC4H10 

      P =   &
          rcemis(NC4H10,k) 
      ! L = 0.0
 

      xnew(NC4H10)=  xold(NC4H10) + dt2 * P

!-> OXYL 

      P =   &
          rcemis(OXYL,k) 
      ! L = 0.0
 

      xnew(OXYL)=  xold(OXYL) + dt2 * P

!-> C5H8 

      P =   &
          rcemis(C5H8,k)  &
        + 0 !Skip bio rate since rcemis exists 
      ! L = 0.0
 

      xnew(C5H8)=  xold(C5H8) + dt2 * P

!-> APINENE 

      P =   &
          rcemis(APINENE,k) 
      ! L = 0.0
 

      xnew(APINENE)=  xold(APINENE) + dt2 * P

!-> C3H6 

      P =   &
          rcemis(C3H6,k) 
      ! L = 0.0
 

      xnew(C3H6)=  xold(C3H6) + dt2 * P

!-> C2H4 

      P =   &
          rcemis(C2H4,k) 
      ! L = 0.0
 

      xnew(C2H4)=  xold(C2H4) + dt2 * P

!-> C2H5OH 

      P =   &
          rcemis(C2H5OH,k) 
      ! L = 0.0
 

      xnew(C2H5OH)=  xold(C2H5OH) + dt2 * P

!-> CH3OH 

      P =   &
          rcemis(CH3OH,k) 
      ! L = 0.0
 

      xnew(CH3OH)=  xold(CH3OH) + dt2 * P

!-> PPM25_FIRE 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for PPM25_FIRE! xnew(PPM25_FIRE)=  max(0.0, xold(PPM25_FIRE))

!-> MEK 

      P =   &
          rcemis(MEK,k) 
      ! L = 0.0
 

      xnew(MEK)=  xold(MEK) + dt2 * P

!-> GLYOX 

      P =   &
          rcemis(GLYOX,k) 
      ! L = 0.0
 

      xnew(GLYOX)=  xold(GLYOX) + dt2 * P

!-> MGLYOX 

      P =   &
          rcemis(MGLYOX,k) 
      ! L = 0.0
 

      xnew(MGLYOX)=  xold(MGLYOX) + dt2 * P

!-> UNREAC 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for UNREAC! xnew(UNREAC)=  max(0.0, xold(UNREAC))

!-> V1702A02B_1 

      P =   &
          rcemis(V1702A02B_1,k) 
      ! L = 0.0
 

      xnew(V1702A02B_1)=  xold(V1702A02B_1) + dt2 * P

!-> V1702A02B_2 

      P =   &
          rcemis(V1702A02B_2,k) 
      ! L = 0.0
 

      xnew(V1702A02B_2)=  xold(V1702A02B_2) + dt2 * P

!-> V1702A02B_3 

      P =   &
          rcemis(V1702A02B_3,k) 
      ! L = 0.0
 

      xnew(V1702A02B_3)=  xold(V1702A02B_3) + dt2 * P

!-> V1702A02B_4 

      P =   &
          rcemis(V1702A02B_4,k) 
      ! L = 0.0
 

      xnew(V1702A02B_4)=  xold(V1702A02B_4) + dt2 * P

!-> V1702A02B_5 

      P =   &
          rcemis(V1702A02B_5,k) 
      ! L = 0.0
 

      xnew(V1702A02B_5)=  xold(V1702A02B_5) + dt2 * P

!-> V1702A02B_6 

      P =   &
          rcemis(V1702A02B_6,k) 
      ! L = 0.0
 

      xnew(V1702A02B_6)=  xold(V1702A02B_6) + dt2 * P

!-> V1702A02B_7 

      P =   &
          rcemis(V1702A02B_7,k) 
      ! L = 0.0
 

      xnew(V1702A02B_7)=  xold(V1702A02B_7) + dt2 * P

!-> V1702A02B_8 

      P =   &
          rcemis(V1702A02B_8,k) 
      ! L = 0.0
 

      xnew(V1702A02B_8)=  xold(V1702A02B_8) + dt2 * P

!-> V1702A02B_9 

      P =   &
          rcemis(V1702A02B_9,k) 
      ! L = 0.0
 

      xnew(V1702A02B_9)=  xold(V1702A02B_9) + dt2 * P

!-> SEASALT_F 

      P =   &
          rcemis(SEASALT_F,k) 
      ! L = 0.0
 

      xnew(SEASALT_F)=  xold(SEASALT_F) + dt2 * P

!-> SEASALT_C 

      P =   &
          rcemis(SEASALT_C,k) 
      ! L = 0.0
 

      xnew(SEASALT_C)=  xold(SEASALT_C) + dt2 * P

!-> DUST_ROAD_F 

      P =   &
          rcemis(DUST_ROAD_F,k) 
      ! L = 0.0
 

      xnew(DUST_ROAD_F)=  xold(DUST_ROAD_F) + dt2 * P

!-> DUST_ROAD_C 

      P =   &
          rcemis(DUST_ROAD_C,k) 
      ! L = 0.0
 

      xnew(DUST_ROAD_C)=  xold(DUST_ROAD_C) + dt2 * P

!-> DUST_WB_F 

      P =   &
          rcemis(DUST_WB_F,k) 
      ! L = 0.0
 

      xnew(DUST_WB_F)=  xold(DUST_WB_F) + dt2 * P

!-> DUST_WB_C 

      P =   &
          rcemis(DUST_WB_C,k) 
      ! L = 0.0
 

      xnew(DUST_WB_C)=  xold(DUST_WB_C) + dt2 * P

!-> DUST_SAH_F 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for DUST_SAH_F! xnew(DUST_SAH_F)=  max(0.0, xold(DUST_SAH_F))

!-> DUST_SAH_C 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for DUST_SAH_C! xnew(DUST_SAH_C)=  max(0.0, xold(DUST_SAH_C))
