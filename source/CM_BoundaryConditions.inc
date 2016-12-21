!---- Boundary conditions -------
!---- SECTION BICs 
  !Available BCs are indexed:
  !use GlobalBCs_ml,  only:  NGLOB_BC  &!  indices from UiO model
  !              ,IBC_SO2, IBC_SO4, IBC_HCHO, IBC_CH3CHO &
  !              ,IBC_O3,IBC_HNO3,IBC_PAN,IBC_CO,IBC_C2H6   &
  !              ,IBC_C4H10, IBC_NO ,IBC_NO2,IBC_NH4_f,IBC_NO3_f,IBC_NO3_c&
  !             ,IBC_H2O2,IBC_CH3COO2

 ! mappings for species from Logan + obs model given with IBC index. 
 ! EMEP model species have IXADV indices.

  bc2xn_adv(IBC_SEASALT_F  ,IXADV_SEASALT_F  )   =   1.0
  bc2xn_adv(IBC_SEASALT_C  ,IXADV_SEASALT_C  )   =   1.0

! Skipped:
!  bc2xn_adv(IBC_SEASALT_G  ,IXADV_SEASALT_G  )   =   1.0
!---- SECTION BICs 
  !Available BCs are indexed:
  !use GlobalBCs_ml,  only:  NGLOB_BC  &!  indices from UiO model
  !              ,IBC_SO2, IBC_SO4, IBC_HCHO, IBC_CH3CHO &
  !              ,IBC_O3,IBC_HNO3,IBC_PAN,IBC_CO,IBC_C2H6   &
  !              ,IBC_C4H10, IBC_NO ,IBC_NO2,IBC_NH4_f,IBC_NO3_f,IBC_NO3_c&
  !             ,IBC_H2O2,IBC_CH3COO2

 ! mappings for species from Logan + obs model given with IBC index. 
 ! EMEP model species have IXADV indices.

  bc2xn_adv(IBC_DUST_f  ,IXADV_DUST_SAH_F )  = 1.0  !Dust
  bc2xn_adv(IBC_DUST_c  ,IXADV_DUST_SAH_C )  = 1.0  !Dust

