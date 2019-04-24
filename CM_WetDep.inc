type(dep_t), public, dimension(NWETDEP_ADV), save :: CM_WDepMap = (/ &
!        species     surrogate      setrate
  dep_t("N2O5      ", "HNO3      ", -999.0 ) &
, dep_t("HNO3      ", "HNO3      ", -999.0 ) &
, dep_t("HONO      ", "HNO3      ", -999.0 ) &
, dep_t("HO2NO2    ", "HNO3      ", -999.0 ) &
, dep_t("HCHO      ", "HCHO      ", -999.0 ) &
, dep_t("H2O2      ", "H2O2      ", -999.0 ) &
, dep_t("SO2       ", "SO2       ", -999.0 ) &
, dep_t("SO4       ", "PMf       ", -999.0 ) &
, dep_t("NH3       ", "NH3       ", -999.0 ) &
, dep_t("NO3_f     ", "PMf       ", -999.0 ) &
, dep_t("NO3_c     ", "PMc       ", -999.0 ) &
, dep_t("NH4_f     ", "PMf       ", -999.0 ) &
, dep_t("SQT_SOA_NV", "PMf       ", -999.0 ) &
, dep_t("POM_f_WOOD", "PMf       ", -999.0 ) &
, dep_t("POM_f_FFUEL", "PMf       ", -999.0 ) &
, dep_t("POM_c_FFUEL", "PMc       ", -999.0 ) &
, dep_t("EC_f_WOOD_new", "ECfn      ", -999.0 ) &
, dep_t("EC_f_WOOD_age", "PMf       ", -999.0 ) &
, dep_t("EC_c_WOOD ", "PMc       ", -999.0 ) &
, dep_t("EC_f_FFUEL_new", "ECfn      ", -999.0 ) &
, dep_t("EC_f_FFUEL_age", "PMf       ", -999.0 ) &
, dep_t("EC_c_FFUEL", "PMc       ", -999.0 ) &
, dep_t("REMPPM25  ", "PMf       ", -999.0 ) &
, dep_t("REMPPM_c  ", "PMc       ", -999.0 ) &
, dep_t("ASOC_ng100", "PMf       ", -999.0 ) &
, dep_t("ASOC_ug1  ", "PMf       ", -999.0 ) &
, dep_t("ASOC_ug10 ", "PMf       ", -999.0 ) &
, dep_t("ASOC_ug1e2", "PMf       ", -999.0 ) &
, dep_t("ASOC_ug1e3", "PMf       ", -999.0 ) &
, dep_t("NON_C_ASOA_ng100", "PMf       ", -999.0 ) &
, dep_t("NON_C_ASOA_ug1", "PMf       ", -999.0 ) &
, dep_t("NON_C_ASOA_ug10", "PMf       ", -999.0 ) &
, dep_t("NON_C_ASOA_ug1e2", "PMf       ", -999.0 ) &
, dep_t("NON_C_ASOA_ug1e3", "PMf       ", -999.0 ) &
, dep_t("BSOC_ng100", "PMf       ", -999.0 ) &
, dep_t("BSOC_ug1  ", "PMf       ", -999.0 ) &
, dep_t("BSOC_ug10 ", "PMf       ", -999.0 ) &
, dep_t("BSOC_ug1e2", "PMf       ", -999.0 ) &
, dep_t("BSOC_ug1e3", "PMf       ", -999.0 ) &
, dep_t("NON_C_BSOA_ng100", "PMf       ", -999.0 ) &
, dep_t("NON_C_BSOA_ug1", "PMf       ", -999.0 ) &
, dep_t("NON_C_BSOA_ug10", "PMf       ", -999.0 ) &
, dep_t("NON_C_BSOA_ug1e2", "PMf       ", -999.0 ) &
, dep_t("NON_C_BSOA_ug1e3", "PMf       ", -999.0 ) &
, dep_t("FFFUEL_ng10", "PMf       ", -999.0 ) &
, dep_t("WOODOA_ng10", "PMf       ", -999.0 ) &
, dep_t("FFIRE_OM  ", "PMf       ", -999.0 ) &
, dep_t("FFIRE_BC  ", "PMf       ", -999.0 ) &
, dep_t("FFIRE_REMPPM25", "PMf       ", -999.0 ) &
, dep_t("SeaSalt_f ", "SSf       ", -999.0 ) &
, dep_t("SeaSalt_c ", "SSc       ", -999.0 ) &
, dep_t("Dust_ROAD_f", "PMf       ", -999.0 ) &
, dep_t("Dust_ROAD_c", "PMc       ", -999.0 ) &
, dep_t("Dust_WB_f ", "PMf       ", -999.0 ) &
, dep_t("Dust_WB_c ", "PMc       ", -999.0 ) &
, dep_t("Dust_SAH_f", "PMf       ", -999.0 ) &
, dep_t("Dust_SAH_c", "PMc       ", -999.0 ) &
, dep_t("Ash_f     ", "PMf       ", -999.0 ) &
, dep_t("Ash_c     ", "PMc       ", -999.0 ) &
/)
