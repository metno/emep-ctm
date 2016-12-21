  integer, parameter, public ::  NEMIS_Specs  = 38
  character(len=11), save, dimension(NEMIS_Specs), public:: &
      EMIS_Specs =  (/ &
             "TRACER1    " &
           , "TRACER2    " &
           , "NO         " &
           , "NO2        " &
           , "SHIPNOX    " &
           , "SO2        " &
           , "SO4        " &
           , "CO         " &
           , "NH3        " &
           , "C2H6       " &
           , "NC4H10     " &
           , "C2H4       " &
           , "C3H6       " &
           , "OXYL       " &
           , "HCHO       " &
           , "CH3CHO     " &
           , "MEK        " &
           , "C2H5OH     " &
           , "CH3OH      " &
           , "GLYOX      " &
           , "MGLYOX     " &
           , "C5H8       " &
           , "APINENE    " &
           , "V1702A02B_1" &
           , "V1702A02B_2" &
           , "V1702A02B_3" &
           , "V1702A02B_4" &
           , "V1702A02B_5" &
           , "V1702A02B_6" &
           , "V1702A02B_7" &
           , "V1702A02B_8" &
           , "V1702A02B_9" &
           , "SEASALT_F  " &
           , "SEASALT_C  " &
           , "DUST_WB_F  " &
           , "DUST_WB_C  " &
           , "DUST_ROAD_F" &
           , "DUST_ROAD_C" &
 /)
 