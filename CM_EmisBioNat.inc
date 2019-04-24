  integer, parameter, public ::  NEMIS_BioNat  = 14
  character(len=12), save, dimension(NEMIS_BioNat), public:: &
      EMIS_BioNat =  (/ &
             "C5H8        " &
           , "BIOTERP     " &
           , "NO          " &
           , "ASH_F       " &
           , "ASH_C       " &
           , "SEASALT_F   " &
           , "SEASALT_C   " &
           , "DUST_WB_F   " &
           , "DUST_WB_C   " &
           , "DUST_ROAD_F " &
           , "DUST_ROAD_C " &
           , "POLLEN_BIRCH" &
           , "POLLEN_OLIVE" &
           , "POLLEN_GRASS" &
 /)
 