  integer, parameter, public ::  NEMIS_BioNat  = 13
  character(len=11), save, dimension(NEMIS_BioNat), public:: &
      EMIS_BioNat =  (/ &
             "C5H8       " &
           , "NO         " &
           , "NH3        " &
           , "ASH_F      " &
           , "ASH_C      " &
           , "BIOTERP    " &
           , "SEASALT_F  " &
           , "SEASALT_C  " &
           , "DUST_WB_F  " &
           , "DUST_WB_C  " &
           , "DUST_ROAD_F" &
           , "DUST_ROAD_C" &
           , "RN222      " &
 /)
 