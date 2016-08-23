integer, public, parameter :: NUM_EXT = 22
type(ExtEffMap), public, dimension(NUM_EXT), parameter :: ExtMap=(/&
  ExtEffMap(DUST_ROAD_C     ,CEXT_DDc ,DRY_MODE),&
  ExtEffMap(DUST_ROAD_F     ,CEXT_DDf ,DRY_MODE),&
  ExtEffMap(DUST_SAH_C      ,CEXT_DDc ,WET_MODE),&
  ExtEffMap(DUST_SAH_F      ,CEXT_DDf ,WET_MODE),&
  ExtEffMap(DUST_WB_C       ,CEXT_DDc ,WET_MODE),&
  ExtEffMap(DUST_WB_F       ,CEXT_DDf ,WET_MODE),&
  ExtEffMap(EC_F_FFUEL_AGE  ,CEXT_ECa ,WET_MODE),&
  ExtEffMap(EC_F_FFUEL_NEW  ,CEXT_ECn ,WET_MODE),&
  ExtEffMap(EC_F_WOOD_AGE   ,CEXT_ECa ,WET_MODE),&
  ExtEffMap(EC_F_WOOD_NEW   ,CEXT_ECn ,WET_MODE),&
  ExtEffMap(FFIRE_BC        ,CEXT_EC  ,WET_MODE),&
  ExtEffMap(FFIRE_OM        ,CEXT_OC  ,WET_MODE),&
  ExtEffMap(FFIRE_REMPPM25  ,CEXT_DDf ,WET_MODE),&
  ExtEffMap(NH4_F           ,CEXT_NH4f,WET_MODE),&
  ExtEffMap(NO3_C           ,CEXT_NO3c,WET_MODE),&
  ExtEffMap(NO3_F           ,CEXT_NO3f,WET_MODE),&
  ExtEffMap(PART_OM_F       ,CEXT_OC  ,WET_MODE),&
  ExtEffMap(REMPPM25        ,CEXT_DDf ,WET_MODE),&
  ExtEffMap(REMPPM_C        ,CEXT_DDc ,WET_MODE),&
  ExtEffMap(SEASALT_C       ,CEXT_SSc ,WET_MODE),&
  ExtEffMap(SEASALT_F       ,CEXT_SSf ,WET_MODE),&
  ExtEffMap(SO4             ,CEXT_SO4 ,WET_MODE)/)
