!-------------------------------------------------------------------------!
!     Mapping of emissions from FINNv1 to EMEP EmChem09 species
!     Included by  ForestFire_ml 
!-------------------------------------------------------------------------!

  character(len=*), public, parameter :: &
      BiomassBurningMapping = "FINNv1toEmChem09"
  integer, private, parameter :: &
    NBB_DEFS  = 19 & ! No mapping lines below
   ,NEMEPSPECS  = 15   ! No EMEP chemical mech specs used

  !-----------------------------------------------------------------------!
  ! Column-1 gives FINNv1 species ,
  ! Column-2 gives MW if needed to get to kg/day.
  !  (Some FINN emissions are given with mass basis, kg/day, so just
  !    set to 1.0, except for OC, which is given in units of C and 
  !    needs to be scaled from OC to OM. Others are given as mole/day, 
  !    so multiply by MW/1000 to get kg/day )
  ! IMPORTANT - use just one value for unitsfac per txt
  ! Column-3 gives mass fraction of FINN emissions assign to EMEP species
  ! Column-4 gives EMEP species 
  !-----------------------------------------------------------------------!

  !     BBname    unitsfac frac  emep:
  type(bbtype), private, dimension(NBB_DEFS) :: FF_defs = (/ &
    bbtype("CO  ",0.028 , 1.000, CO          ) &
   ,bbtype("NO  ",0.030 , 1.000, NO          ) &
   ,bbtype("NO2 ",0.046 , 1.000, NO2         ) &
   ,bbtype("SO2 ",0.064 , 1.000, SO2         ) &
   ,bbtype("NH3 ",0.017 , 1.000, NH3         ) &
   ,bbtype("ACET",0.058 , 1.000, C2H6        ) & ! acetone
   ,bbtype("ALD2",0.044 , 1.000, CH3CHO      ) &
   ,bbtype("ALK4",0.058 , 1.000, NC4H10      ) &
   ,bbtype("C2H6",0.030 , 1.000, C2H6        ) &
   ,bbtype("C3H8",0.044 , 0.700, NC4H10      ) & ! obs
   ,bbtype("CH2O",0.030 , 1.000, HCHO        ) &
   ,bbtype("ISOP",0.068 , 1.000, C5H8        ) &
   ,bbtype("MEK ",0.072 , 1.000, MEK         ) &
   ,bbtype("PRPE",0.042 , 1.000, C3H6        ) & !CHECK!
! We read in OC and PM25, but want OM and REMPPM25
   ,bbtype("PM25",1.0   , 1.000, FFIRE_REMPPM25 ) & !  Will need to subtract OM, BC
   ,bbtype("OC  ",1.7   , 1.000, FFIRE_OM    ) & ! Put OM/OC=1.7 in fac
   ,bbtype("BC  ",1.0   , 1.000, FFIRE_BC    ) &
! Subtract, assuming OM/OC=1.7. ForestFire_ml will pevent zeros
   ,bbtype("OC  ",-1.7  , 1.000, FFIRE_REMPPM25 ) & ! Will subtract OM
   ,bbtype("BC  ",-1.0  , 1.000, FFIRE_REMPPM25 ) & ! Will subtract BC
 /)
! Not used in EmChem09 standard
 !,bbtype("CH4  "      , 1.000, CH4         ) &
 ! ,bbtype("HCN "       , 1.000, HCN         ) &
