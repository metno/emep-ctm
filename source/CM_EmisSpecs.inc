  integer, parameter ::  NEMIS_SPECS  = 22
  character(len=12), save, dimension(NEMIS_SPECS):: &
      EMIS_SPECS =  (/ &
             "NO          " &
           , "NO2         " &
           , "SO2         " &
           , "SO4         " &
           , "NH3         " &
           , "CO          " &
           , "C2H6        " &
           , "NC4H10      " &
           , "C2H4        " &
           , "C3H6        " &
           , "OXYL        " &
           , "HCHO        " &
           , "CH3CHO      " &
           , "MEK         " &
           , "C2H5OH      " &
           , "CH3OH       " &
           , "GLYOX       " &
           , "MGLYOX      " &
           , "C5H8        " &
           , "PPM25       " &
           , "PPM_C       " &
           , "PPM25_FIRE  " &
 /)
 