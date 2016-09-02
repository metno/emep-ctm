integer, public, parameter :: NUM_EXT = 24
type(ExtEffMap), public, dimension(NUM_EXT), parameter :: ExtMap=(/&
  ExtEffMap(ASH_C           ,CEXT_DDc ),&
  ExtEffMap(ASH_F           ,CEXT_DDf ),&
  ExtEffMap(DUST_ROAD_C     ,CEXT_DDc ),&
  ExtEffMap(DUST_ROAD_F     ,CEXT_DDf ),&
  ExtEffMap(DUST_SAH_C      ,CEXT_DDc ),&
  ExtEffMap(DUST_SAH_F      ,CEXT_DDf ),&
  ExtEffMap(DUST_WB_C       ,CEXT_DDc ),&
  ExtEffMap(DUST_WB_F       ,CEXT_DDf ),&
  ExtEffMap(EC_F_FFUEL_AGE  ,CEXT_ECa ),&
  ExtEffMap(EC_F_FFUEL_NEW  ,CEXT_ECn ),&
  ExtEffMap(EC_F_WOOD_AGE   ,CEXT_ECa ),&
  ExtEffMap(EC_F_WOOD_NEW   ,CEXT_ECn ),&
  ExtEffMap(FFIRE_BC        ,CEXT_EC  ),&
  ExtEffMap(FFIRE_OM        ,CEXT_OC  ),&
  ExtEffMap(FFIRE_REMPPM25  ,CEXT_DDf ),&
  ExtEffMap(NH4_F           ,CEXT_NH4f),&
  ExtEffMap(NO3_C           ,CEXT_NO3c),&
  ExtEffMap(NO3_F           ,CEXT_NO3f),&
  ExtEffMap(OM25_P          ,CEXT_OC  ),&
  ExtEffMap(REMPPM25        ,CEXT_DDf ),&
  ExtEffMap(REMPPM_C        ,CEXT_DDc ),&
  ExtEffMap(SEASALT_C       ,CEXT_SSc ),&
  ExtEffMap(SEASALT_F       ,CEXT_SSf ),&
  ExtEffMap(SO4             ,CEXT_SO4 )/)
