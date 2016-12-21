  integer, public, parameter ::  NDRYDEP_ADV  = 31
  type(depmap), public, dimension(NDRYDEP_ADV), parameter:: DDepMap= (/ &
        depmap( IXADV_O3, CDDEP_O3, -1) & 
      , depmap( IXADV_NO2, CDDEP_NO2, -1) & 
      , depmap( IXADV_SHIPNOX, CDDEP_NO2, -1) & 
      , depmap( IXADV_HNO3, CDDEP_HNO3, -1) & 
      , depmap( IXADV_CH3CHO, CDDEP_ALD, -1) & 
      , depmap( IXADV_PAN, CDDEP_PAN, -1) & 
      , depmap( IXADV_SO2, CDDEP_SO2, -1) & 
      , depmap( IXADV_SO4, CDDEP_PMfS, -1) & 
      , depmap( IXADV_NH3, CDDEP_NH3, -1) & 
      , depmap( IXADV_NO3_f, CDDEP_PMfN, -1) & 
      , depmap( IXADV_NO3_c, CDDEP_PMc, -1) & 
      , depmap( IXADV_NH4_f, CDDEP_PMfN, -1) & 
      , depmap( IXADV_HCHO, CDDEP_HCHO, -1) & 
      , depmap( IXADV_PPM25_FIRE, CDDEP_PMfS, -1) & 
      , depmap( IXADV_V1702A02B_1, CDDEP_PMc, -1) & 
      , depmap( IXADV_V1702A02B_2, CDDEP_PMc, -1) & 
      , depmap( IXADV_V1702A02B_3, CDDEP_PMc, -1) & 
      , depmap( IXADV_V1702A02B_4, CDDEP_PMc, -1) & 
      , depmap( IXADV_V1702A02B_5, CDDEP_PMc, -1) & 
      , depmap( IXADV_V1702A02B_6, CDDEP_PMc, -1) & 
      , depmap( IXADV_V1702A02B_7, CDDEP_PMc, -1) & 
      , depmap( IXADV_V1702A02B_8, CDDEP_PMc, -1) & 
      , depmap( IXADV_V1702A02B_9, CDDEP_PMc, -1) & 
      , depmap( IXADV_SeaSalt_f, CDDEP_PMfS, -1) & 
      , depmap( IXADV_SeaSalt_c, CDDEP_SSc, -1) & 
      , depmap( IXADV_DUST_ROAD_F, CDDEP_PMfS, -1) & 
      , depmap( IXADV_DUST_ROAD_C, CDDEP_DUc, -1) & 
      , depmap( IXADV_Dust_wb_f, CDDEP_PMfS, -1) & 
      , depmap( IXADV_Dust_wb_c, CDDEP_DUc, -1) & 
      , depmap( IXADV_Dust_sah_f, CDDEP_PMfS, -1) & 
      , depmap( IXADV_Dust_sah_c, CDDEP_DUc, -1) & 
  /)
