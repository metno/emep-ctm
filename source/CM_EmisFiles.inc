  integer, parameter, public ::  NEMIS_File  = 2
  character(len=3), save, dimension(NEMIS_File), public:: &
      EMIS_File =  (/ &
             "nox" &
           , "voc" &
 /)
 