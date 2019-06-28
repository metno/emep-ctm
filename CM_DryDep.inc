type(dep_t), public, dimension(NDRYDEP_ADV), save :: CM_DDepMap = (/ &
!        species     surrogate      setrate
  dep_t("O3        ", "O3        ", -999.0 ) &
, dep_t("NO2       ", "NO2       ", -999.0 ) &
, dep_t("N2O5      ", "N2O5      ", -999.0 ) &
, dep_t("H2O2      ", "H2O2      ", -999.0 ) &
, dep_t("HONO      ", "HNO2      ", -999.0 ) &
, dep_t("HNO3      ", "HNO3      ", -999.0 ) &
, dep_t("HO2NO2    ", "HO2NO2    ", -999.0 ) &
, dep_t("HCHO      ", "HCHO      ", -999.0 ) &
, dep_t("CH3CHO    ", "ALD       ", -999.0 ) &
, dep_t("MACR      ", "MEK       ", -999.0 ) &
, dep_t("MEK       ", "MEK       ", -999.0 ) &
, dep_t("ACETOL    ", "ACETOL    ", -999.0 ) &
, dep_t("GLYOX     ", "GLYOX     ", -999.0 ) &
, dep_t("MGLYOX    ", "MGLYOX    ", -999.0 ) &
, dep_t("BIACET    ", "MEK       ", -999.0 ) &
, dep_t("C5DICARB  ", "C5DICARB  ", -999.0 ) &
, dep_t("CH3OOH    ", "MEOOH     ", -999.0 ) &
, dep_t("C2H5OOH   ", "MEOOH     ", -999.0 ) &
, dep_t("BURO2H    ", "C3H7OOH   ", -999.0 ) &
, dep_t("ETRO2H    ", "SHISOLOOH ", -999.0 ) &
, dep_t("PRRO2H    ", "SHISOLOOH ", -999.0 ) &
, dep_t("MEKO2H    ", "MDSOLOOH  ", -999.0 ) &
, dep_t("ISRO2H    ", "SHISOLOOH ", -999.0 ) &
, dep_t("C5DICAROOH", "VHISOLNO3 ", -999.0 ) &
, dep_t("HPALD     ", "MDSOLOOH  ", -999.0 ) &
, dep_t("MACROOH   ", "VHISOLOOH ", -999.0 ) &
, dep_t("OXYO2H    ", "VHISOLOOH ", -999.0 ) &
, dep_t("CH3CO3H   ", "PAA       ", -999.0 ) &
, dep_t("PACALD    ", "CO2C5OH   ", -999.0 ) &
, dep_t("IEPOX     ", "HISOLF0   ", -999.0 ) &
, dep_t("SC4H9NO3  ", "VLSOLNO3  ", -999.0 ) &
, dep_t("NALD      ", "PAN       ", -999.0 ) &
, dep_t("ISON      ", "MDNO3OH   ", -999.0 ) &
, dep_t("PAN       ", "PAN       ", -999.0 ) &
, dep_t("MPAN      ", "PAN       ", -999.0 ) &
, dep_t("SO2       ", "SO2       ", -999.0 ) &
, dep_t("SO4       ", "PMf       ", -999.0 ) &
, dep_t("NH3       ", "NH3       ", -999.0 ) &
, dep_t("NO3_f     ", "PMfNO3    ", -999.0 ) &
, dep_t("NO3_c     ", "PMc       ", -999.0 ) &
, dep_t("NH4_f     ", "PMfNH4    ", -999.0 ) &
, dep_t("POM_f_wood", "PMf       ", -999.0 ) &
, dep_t("POM_c_wood", "PMc       ", -999.0 ) &
, dep_t("POM_f_ffuel", "PMf       ", -999.0 ) &
, dep_t("POM_c_ffuel", "PMc       ", -999.0 ) &
, dep_t("EC_f_wood_new", "PMf       ", -999.0 ) &
, dep_t("EC_f_wood_age", "PMf       ", -999.0 ) &
, dep_t("EC_c_wood ", "PMc       ", -999.0 ) &
, dep_t("EC_f_ffuel_new", "PMf       ", -999.0 ) &
, dep_t("EC_f_ffuel_age", "PMf       ", -999.0 ) &
, dep_t("EC_c_ffuel", "PMc       ", -999.0 ) &
, dep_t("pSO4f     ", "PMf       ", -999.0 ) &
, dep_t("pSO4c     ", "PMc       ", -999.0 ) &
, dep_t("remPPM25  ", "PMf       ", -999.0 ) &
, dep_t("remPPM_c  ", "PMc       ", -999.0 ) &
, dep_t("ASOC_ng1e2", "ALD       ", -999.0 ) &
, dep_t("ASOC_ug1  ", "ALD       ", -999.0 ) &
, dep_t("ASOC_ug10 ", "ALD       ", -999.0 ) &
, dep_t("ASOC_ug1e2", "ALD       ", -999.0 ) &
, dep_t("ASOC_ug1e3", "ALD       ", -999.0 ) &
, dep_t("non_C_ASOA_ng1e2", "ALD       ", -999.0 ) &
, dep_t("non_C_ASOA_ug1", "ALD       ", -999.0 ) &
, dep_t("non_C_ASOA_ug10", "ALD       ", -999.0 ) &
, dep_t("non_C_ASOA_ug1e2", "ALD       ", -999.0 ) &
, dep_t("non_C_ASOA_ug1e3", "ALD       ", -999.0 ) &
, dep_t("BSOC_ng1e2", "ALD       ", -999.0 ) &
, dep_t("BSOC_ug1  ", "ALD       ", -999.0 ) &
, dep_t("BSOC_ug10 ", "ALD       ", -999.0 ) &
, dep_t("BSOC_ug1e2", "ALD       ", -999.0 ) &
, dep_t("BSOC_ug1e3", "ALD       ", -999.0 ) &
, dep_t("non_C_BSOA_ng1e2", "ALD       ", -999.0 ) &
, dep_t("non_C_BSOA_ug1", "ALD       ", -999.0 ) &
, dep_t("non_C_BSOA_ug10", "ALD       ", -999.0 ) &
, dep_t("non_C_BSOA_ug1e2", "ALD       ", -999.0 ) &
, dep_t("non_C_BSOA_ug1e3", "ALD       ", -999.0 ) &
, dep_t("SQT_SOA_nv", "PMf       ", -999.0 ) &
, dep_t("TERPOOH   ", "LHISOLOOH ", -999.0 ) &
, dep_t("MVK       ", "MEK       ", -999.0 ) &
, dep_t("Ash_f     ", "PMf       ", -999.0 ) &
, dep_t("Ash_c     ", "PMc       ", -999.0 ) &
, dep_t("shipNOx   ", "NO2       ", -999.0 ) &
, dep_t("ffire_OM  ", "PMf       ", -999.0 ) &
, dep_t("ffire_BC  ", "PMf       ", -999.0 ) &
, dep_t("ffire_remPPM25", "PMf       ", -999.0 ) &
, dep_t("SeaSalt_f ", "SSf       ", -999.0 ) &
, dep_t("SeaSalt_c ", "SSc       ", -999.0 ) &
, dep_t("Dust_road_f", "DUf       ", -999.0 ) &
, dep_t("Dust_road_c", "DUc       ", -999.0 ) &
, dep_t("Dust_wb_f ", "DUf       ", -999.0 ) &
, dep_t("Dust_wb_c ", "DUc       ", -999.0 ) &
, dep_t("Dust_sah_f", "DUf       ", -999.0 ) &
, dep_t("Dust_sah_c", "DUc       ", -999.0 ) &
/)
