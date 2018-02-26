  integer, parameter, public ::  NEMIS_Specs  = 44
  character(len=14), save, dimension(NEMIS_Specs), public:: &
      EMIS_Specs =  (/ &
             "NO            " &
           , "NO2           " &
           , "SHIPNOX       " &
           , "SO2           " &
           , "SO4           " &
           , "CO            " &
           , "NH3           " &
           , "C2H6          " &
           , "NC4H10        " &
           , "C2H4          " &
           , "C3H6          " &
           , "OXYL          " &
           , "HCHO          " &
           , "CH3CHO        " &
           , "MEK           " &
           , "C2H5OH        " &
           , "CH3OH         " &
           , "GLYOX         " &
           , "MGLYOX        " &
           , "C5H8          " &
           , "ASH_F         " &
           , "ASH_C         " &
           , "FFIRE_CO      " &
           , "BIOTERP       " &
           , "POM_F_FFUEL   " &
           , "POM_C_FFUEL   " &
           , "EC_F_FFUEL_NEW" &
           , "EC_F_FFUEL_AGE" &
           , "EC_C_FFUEL    " &
           , "POM_F_WOOD    " &
           , "EC_F_WOOD_NEW " &
           , "EC_F_WOOD_AGE " &
           , "REMPPM25      " &
           , "REMPPM_C      " &
           , "FFIRE_OM      " &
           , "FFIRE_BC      " &
           , "FFIRE_REMPPM25" &
           , "SEASALT_F     " &
           , "SEASALT_C     " &
           , "DUST_WB_F     " &
           , "DUST_WB_C     " &
           , "DUST_ROAD_F   " &
           , "DUST_ROAD_C   " &
           , "RN222         " &
 /)
 