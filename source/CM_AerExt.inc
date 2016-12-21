integer, public, parameter :: NUM_EXT = 22
type(ExtEffMap), public, dimension(NUM_EXT), parameter :: ExtMap=(/&
  ExtEffMap(DUST_ROAD_C     ,CEXT_DDc ),&
  ExtEffMap(DUST_ROAD_F     ,CEXT_DDf ),&
  ExtEffMap(DUST_SAH_C      ,CEXT_DDc ),&
  ExtEffMap(DUST_SAH_F      ,CEXT_DDf ),&
  ExtEffMap(DUST_WB_C       ,CEXT_DDc ),&
  ExtEffMap(DUST_WB_F       ,CEXT_DDf ),&
  ExtEffMap(NH4_F           ,CEXT_NH4f),&
  ExtEffMap(NO3_C           ,CEXT_NO3c),&
  ExtEffMap(NO3_F           ,CEXT_NO3f),&
  ExtEffMap(PPM25_FIRE      ,CEXT_DDf ),&
  ExtEffMap(SEASALT_C       ,CEXT_SSc ),&
  ExtEffMap(SEASALT_F       ,CEXT_SSf ),&
  ExtEffMap(SO4             ,CEXT_SO4 ),&
  ExtEffMap(V1702A02B_1     ,CEXT_DDf ),&
  ExtEffMap(V1702A02B_2     ,CEXT_DDf ),&
  ExtEffMap(V1702A02B_3     ,CEXT_DDf ),&
  ExtEffMap(V1702A02B_4     ,CEXT_DDf ),&
  ExtEffMap(V1702A02B_5     ,CEXT_DDc ),&
  ExtEffMap(V1702A02B_6     ,CEXT_DDc ),&
  ExtEffMap(V1702A02B_7     ,CEXT_DDc ),&
  ExtEffMap(V1702A02B_8     ,CEXT_DDc ),&
  ExtEffMap(V1702A02B_9     ,CEXT_DDc )/)
