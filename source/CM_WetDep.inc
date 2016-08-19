  integer, public, parameter ::  NWETDEP_ADV  = 19
  type(depmap), public, dimension(NWETDEP_ADV), parameter:: WDepMap= (/ &
       depmap( HNO3, CWDEP_HNO3, -1) & 
      , depmap( HONO, CWDEP_HNO3, -1) & 
      , depmap( HCHO, CWDEP_HCHO, -1) & 
      , depmap( H2O2, CWDEP_H2O2, -1) & 
      , depmap( SO2, CWDEP_SO2, -1) & 
      , depmap( SO4, CWDEP_SO4, -1) & 
      , depmap( NH3, CWDEP_NH3, -1) & 
      , depmap( NO3_f, CWDEP_PMf, -1) & 
      , depmap( NO3_c, CWDEP_PMc, -1) & 
      , depmap( NH4_f, CWDEP_PMf, -1) & 
      , depmap( PPM25, CWDEP_PMf, -1) & 
      , depmap( PPM25_FIRE, CWDEP_PMf, -1) & 
      , depmap( PPM_c, CWDEP_PMc, -1) & 
      , depmap( SeaSalt_f, CWDEP_SSf, -1) & 
      , depmap( SeaSalt_c, CWDEP_SSc, -1) & 
      , depmap( SeaSalt_g, CWDEP_SSg, -1) & 
      , depmap( Dust_nat_f, CWDEP_PMf, -1) & 
      , depmap( Dust_nat_c, CWDEP_PMc, -1) & 
      , depmap( Pb210, CWDEP_PMf, -1) & 
  /)
