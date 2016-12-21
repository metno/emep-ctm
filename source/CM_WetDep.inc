  integer, public, parameter ::  NWETDEP_ADV  = 26
  type(depmap), public, dimension(NWETDEP_ADV), parameter:: WDepMap= (/ &
        depmap( IXADV_HNO3, CWDEP_HNO3, -1) & 
      , depmap( IXADV_SO2, CWDEP_SO2, -1) & 
      , depmap( IXADV_SO4, CWDEP_SO4, -1) & 
      , depmap( IXADV_NH3, CWDEP_NH3, -1) & 
      , depmap( IXADV_NO3_f, CWDEP_PMf, -1) & 
      , depmap( IXADV_NO3_c, CWDEP_PMc, -1) & 
      , depmap( IXADV_NH4_f, CWDEP_PMf, -1) & 
      , depmap( IXADV_HCHO, CWDEP_HCHO, -1) & 
      , depmap( IXADV_PPM25_FIRE, CWDEP_PMf, -1) & 
      , depmap( IXADV_V1702A02B_1, CWDEP_PMf, -1) & 
      , depmap( IXADV_V1702A02B_2, CWDEP_PMf, -1) & 
      , depmap( IXADV_V1702A02B_3, CWDEP_PMf, -1) & 
      , depmap( IXADV_V1702A02B_4, CWDEP_PMf, -1) & 
      , depmap( IXADV_V1702A02B_5, CWDEP_PMc, -1) & 
      , depmap( IXADV_V1702A02B_6, CWDEP_PMc, -1) & 
      , depmap( IXADV_V1702A02B_7, CWDEP_PMc, -1) & 
      , depmap( IXADV_V1702A02B_8, CWDEP_PMc, -1) & 
      , depmap( IXADV_V1702A02B_9, CWDEP_PMc, -1) & 
      , depmap( IXADV_SeaSalt_f, CWDEP_SSf, -1) & 
      , depmap( IXADV_SeaSalt_c, CWDEP_SSc, -1) & 
      , depmap( IXADV_DUST_ROAD_F, CWDEP_PMf, -1) & 
      , depmap( IXADV_DUST_ROAD_C, CWDEP_PMc, -1) & 
      , depmap( IXADV_Dust_wb_f, CWDEP_PMf, -1) & 
      , depmap( IXADV_Dust_wb_c, CWDEP_PMc, -1) & 
      , depmap( IXADV_Dust_sah_f, CWDEP_PMf, -1) & 
      , depmap( IXADV_Dust_sah_c, CWDEP_PMc, -1) & 
  /)
