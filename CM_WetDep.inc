type(dep_t), public, dimension(NWETDEP_ADV), save :: CM_WDepMap = (/ &
!        species     surrogate      setrate
  dep_t("N2O5      ", "HNO3      ", -999.0 ) &
, dep_t("H2O2      ", "H2O2      ", -999.0 ) &
, dep_t("HONO      ", "H2O2      ", -999.0 ) &
, dep_t("HNO3      ", "HNO3      ", -999.0 ) &
, dep_t("HO2NO2    ", "0p8       ", -999.0 ) &
, dep_t("HCHO      ", "ROOH      ", -999.0 ) &
, dep_t("ACETOL    ", "HCHO      ", -999.0 ) &
, dep_t("GLYOX     ", "H2O2      ", -999.0 ) &
, dep_t("MGLYOX    ", "0p3       ", -999.0 ) &
, dep_t("C5DICARB  ", "H2O2      ", -999.0 ) &
, dep_t("ETRO2H    ", "H2O2      ", -999.0 ) &
, dep_t("MEKO2H    ", "1p3       ", -999.0 ) &
, dep_t("ISRO2H    ", "H2O2      ", -999.0 ) &
, dep_t("C5DICAROOH", "H2O2      ", -999.0 ) &
, dep_t("HPALD     ", "H2O2      ", -999.0 ) &
, dep_t("MACROOH   ", "H2O2      ", -999.0 ) &
, dep_t("OXYO2H    ", "H2O2      ", -999.0 ) &
, dep_t("CH3CO3H   ", "ROOH      ", -999.0 ) &
, dep_t("PACALD    ", "H2O2      ", -999.0 ) &
, dep_t("IEPOX     ", "HNO3      ", -999.0 ) &
, dep_t("NALD      ", "ROOH      ", -999.0 ) &
, dep_t("ISON      ", "0p7       ", -999.0 ) &
, dep_t("TERPOOH   ", "H2O2      ", -999.0 ) &
, dep_t("SO2       ", "SO2       ", -999.0 ) &
, dep_t("SO4       ", "PMf       ", -999.0 ) &
, dep_t("NH3       ", "NH3       ", -999.0 ) &
, dep_t("NO3_f     ", "PMf       ", -999.0 ) &
, dep_t("NO3_c     ", "PMc       ", -999.0 ) &
, dep_t("NH4_f     ", "PMf       ", -999.0 ) &
, dep_t("ASOC_ng1e2", "ROOH      ", -999.0 ) &
, dep_t("ASOC_ug1  ", "ROOH      ", -999.0 ) &
, dep_t("ASOC_ug10 ", "ROOH      ", -999.0 ) &
, dep_t("ASOC_ug1e2", "ROOH      ", -999.0 ) &
, dep_t("ASOC_ug1e3", "ROOH      ", -999.0 ) &
, dep_t("non_C_ASOA_ng1e2", "ROOH      ", -999.0 ) &
, dep_t("non_C_ASOA_ug1", "ROOH      ", -999.0 ) &
, dep_t("non_C_ASOA_ug10", "ROOH      ", -999.0 ) &
, dep_t("non_C_ASOA_ug1e2", "ROOH      ", -999.0 ) &
, dep_t("non_C_ASOA_ug1e3", "ROOH      ", -999.0 ) &
, dep_t("BSOC_ng1e2", "ROOH      ", -999.0 ) &
, dep_t("BSOC_ug1  ", "ROOH      ", -999.0 ) &
, dep_t("BSOC_ug10 ", "ROOH      ", -999.0 ) &
, dep_t("BSOC_ug1e2", "ROOH      ", -999.0 ) &
, dep_t("BSOC_ug1e3", "ROOH      ", -999.0 ) &
, dep_t("non_C_BSOA_ng1e2", "ROOH      ", -999.0 ) &
, dep_t("non_C_BSOA_ug1", "ROOH      ", -999.0 ) &
, dep_t("non_C_BSOA_ug10", "ROOH      ", -999.0 ) &
, dep_t("non_C_BSOA_ug1e2", "ROOH      ", -999.0 ) &
, dep_t("non_C_BSOA_ug1e3", "ROOH      ", -999.0 ) &
, dep_t("ffire_OM  ", "PMf       ", -999.0 ) &
, dep_t("ffire_BC  ", "PMf       ", -999.0 ) &
, dep_t("ffire_remPPM25", "PMf       ", -999.0 ) &
, dep_t("ffire_c   ", "PMc       ", -999.0 ) &
, dep_t("SeaSalt_f ", "SSf       ", -999.0 ) &
, dep_t("SeaSalt_c ", "SSc       ", -999.0 ) &
, dep_t("Dust_road_f", "PMf       ", -999.0 ) &
, dep_t("Dust_road_c", "PMc       ", -999.0 ) &
, dep_t("Dust_wb_f ", "PMf       ", -999.0 ) &
, dep_t("Dust_wb_c ", "PMc       ", -999.0 ) &
, dep_t("Dust_sah_f", "PMf       ", -999.0 ) &
, dep_t("Dust_sah_c", "PMc       ", -999.0 ) &
, dep_t("Ash_f     ", "PMf       ", -999.0 ) &
, dep_t("Ash_c     ", "PMc       ", -999.0 ) &
, dep_t("POM_f_Res ", "PMf       ", -999.0 ) &
, dep_t("POM_c_Res ", "PMc       ", -999.0 ) &
, dep_t("POM_f_nonRes", "PMf       ", -999.0 ) &
, dep_t("POM_c_nonRes", "PMc       ", -999.0 ) &
, dep_t("EC_f_Res_new", "ECfn      ", -999.0 ) &
, dep_t("EC_f_Res_age", "PMf       ", -999.0 ) &
, dep_t("EC_c_Res  ", "PMc       ", -999.0 ) &
, dep_t("EC_f_nonRes_new", "ECfn      ", -999.0 ) &
, dep_t("EC_f_nonRes_age", "PMf       ", -999.0 ) &
, dep_t("EC_c_nonRes", "PMc       ", -999.0 ) &
, dep_t("remPPM25_nonRes", "PMf       ", -999.0 ) &
, dep_t("remPPM25_Res", "PMf       ", -999.0 ) &
, dep_t("remPPM_c_nonRes", "PMc       ", -999.0 ) &
, dep_t("remPPM_c_Res", "PMc       ", -999.0 ) &
, dep_t("SQT_SOA_NV", "PMf       ", -999.0 ) &
, dep_t("POLLEN_BIRCH", "POLLw     ", -999.0 ) &
, dep_t("POLLEN_OLIVE", "POLLw     ", -999.0 ) &
, dep_t("POLLEN_RWEED", "POLLw     ", -999.0 ) &
, dep_t("POLLEN_GRASS", "POLLw     ", -999.0 ) &
/)