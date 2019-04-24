  !Available BCs are indexed:
  !use GlobalBCs_ml,  only:  NGLOB_BC  &!  indices from UiO model
  !              ,IBC_SO2, IBC_SO4, IBC_HCHO, IBC_CH3CHO &
  !              ,IBC_O3,IBC_HNO3,IBC_PAN,IBC_CO,IBC_C2H6   &
  !              ,IBC_C4H10, IBC_NO ,IBC_NO2,IBC_NH4_f,IBC_NO3_f,IBC_NO3_c&
  !             ,IBC_H2O2,IBC_CH3COO2

 ! mappings for species from Logan + obs model given with IBC index. 
 ! EMEP model species have IXADV indices.
 ! Changes for CB05 indicated with !CB05

  bc2xn_adv(IBC_O3      ,IXADV_O3      )   =   1.0
  bc2xn_adv(IBC_HNO3    ,IXADV_HNO3    )   =   1.0
  bc2xn_adv(IBC_SO2     ,IXADV_SO2     )   =   1.0
  bc2xn_adv(IBC_SO4     ,IXADV_SO4    )    =   1.0
  bc2xn_adv(IBC_PAN     ,IXADV_PAN     )   =   1.0
  bc2xn_adv(IBC_CO      ,IXADV_CO      )   =   1.0
!BIC  bc2xn_adv(IBC_CO      ,IXADV_CO_BIC  )   =   1.0
!EMEP:
!BIC  bc2xn_adv(IBC_C2H6    ,IXADV_C2H6_BIC)   =   1.0
  bc2xn_adv(IBC_C2H6    ,IXADV_C2H6)       =   1.0
!CB  bc2xn_adv(IBC_C2H6    ,IXADV_ETHA    )   =   1.0 !CB05  C2H6
!EMEP:
  bc2xn_adv(IBC_C4H10   ,IXADV_NC4H10  )   =   1.0
!CB  bc2xn_adv(IBC_C4H10   ,IXADV_PAR     )   =   1.0 !CB05 or 3.0 ????
  bc2xn_adv(IBC_NO      ,IXADV_NO      )   =   1.0
  bc2xn_adv(IBC_NO2     ,IXADV_NO2     )   =   1.0
  bc2xn_adv(IBC_HCHO    ,IXADV_HCHO    )   =   1.0
!EMEP:
!FEB2012  bc2xn_adv(IBC_SEASALT_F  ,IXADV_SEASALT_F  )   =   1.0
!FEB2012  bc2xn_adv(IBC_SEASALT_C  ,IXADV_SEASALT_C  )   =   1.0
!May26 re-instated!
  bc2xn_adv(IBC_SEASALT_F  ,IXADV_SEASALT_F  )   =   1.0
  bc2xn_adv(IBC_SEASALT_C  ,IXADV_SEASALT_C  )   =   1.0
!  bc2xn_adv(IBC_SEASALT_G  ,IXADV_SEASALT_G  )   =   1.0
!CB: bc2xn_adv(IBC_CH3CHO  ,IXADV_ALD2    )   =   1.0 !CB05 ????
  bc2xn_adv(IBC_CH3CHO   ,IXADV_CH3CHO  )   =   1.0
  bc2xn_adv(IBC_NO3_f    ,IXADV_NO3_f  )   =   1.0
  bc2xn_adv(IBC_NO3_c    ,IXADV_NO3_c  )   =   1.0
  bc2xn_adv(IBC_NH4_f    ,IXADV_NH4_f    )   =   1.0
  bc2xn_adv(IBC_H2O2    ,IXADV_H2O2    )   =   1.0
!EMEP:
!SS  bc2xn_adv(IBC_CH3COO2 ,IXADV_CH3COO2 )   =   1.0
!FEB2012  bc2xn_adv(IBC_DUST_f  ,IXADV_DUST_SAH_F )  = 1.0  !Dust
!FEB2012  bc2xn_adv(IBC_DUST_c  ,IXADV_DUST_SAH_C )  = 1.0  !Dust
 bc2xn_adv(IBC_DUST_f  ,IXADV_DUST_SAH_F )  = 1.0  !Dust
 bc2xn_adv(IBC_DUST_c  ,IXADV_DUST_SAH_C )  = 1.0  !Dust
!SKIP  bc2xn_adv(IBC_CH3COO2 ,IXADV_C2O3 )   =   1.0 !CB05
!QUERY??


! The following species are excluded either because they have no corresponding
! species in the emep model, or because they have lifetimes which are so
! short that initialisation is uncessary.
!-----------------------------------------------------------------------------
!u3  bc2xn_adv(IBC_C2H4    ,IXADV_C2H4  )   =   1.0
!u3  bc2xn_adv(IBC_C3H6    ,IXADV_NC4H10  )   =   0.75
!u3  bc2xn_adv(IBC_C6H14   ,IXADV_NC4H10  )   =   1.5   ! ds - scale by C6/C4
!u3  bc2xn_adv(IBC_CH3O2H  ,IXADV_CH3O2H  )   =   1.0
!u3  bc2xn_adv(IBC_ISOPRENE,IXADV_ISOP    )   =   1.0   ! ds-rename
!u3  bc2xn_adv(IBC_RCOHCO  ,IXADV_CH3CHO  )   =   1.0   ! ds - unknown scale
!u3  bc2xn_adv(IBC_CH4     ,IXADV_CH4     )   =   1.0   ! Re-included for DSMACH
!u3  bc2xn_adv(IBC_C3H8    ,IXADV_NC4H10 )   =   0.75   ! mini
!u3  bc2xn_adv(IBC_C3H8    ,IXADV_NC4H10  )   =   0.5   ! ds-split C2H6/NC4H10
!!bc2xn_adv(IBC_NOX ,IXADV_NOX     )   =  -1.0   ! Excluded, we have NO and NO2
!!bc2xn_adv(IBC_C6HXR   ,IXADV_C6HXR   )   =   1.0
!!bc2xn_adv(IBC_HO2NO2  ,IXADV_HO2NO2  )   =   1.0
!!bc2xn_adv(IBC_CH3COY  ,IXADV_CH3COY  )   =   1.0
!!bc2xn_adv(IBC_CH3COX  ,IXADV_CH3COX  )   =   1.0
!!bc2xn_adv(IBC_HO2     ,IXADV_HO2     )   =  -1.0   ! SHort-lived
!!bc2xn_adv(IBC_CH2O2OH ,IXADV_CH2O2OH )   =   1.0
!!bc2xn_adv(IBC_CH3COB  ,IXADV_CH3COB  )   =  -1.0   ! ???
!!bc2xn_adv(IBC_CH3XX   ,IXADV_CH3XX   )   =   1.0
!!bc2xn_adv(IBC_AR1     ,IXADV_AR1     )   =  -1.0   ! ???
!!bc2xn_adv(IBC_AR2     ,IXADV_AR2     )   =   1.0
!!bc2xn_adv(IBC_AR3     ,IXADV_AR3     )   =  -1.0   ! ???
!!bc2xn_adv(IBC_ISOR1   ,IXADV_ISOR1   )   =  -1.0   ! Re-Excluded
!!bc2xn_adv(IBC_ISOK    ,IXADV_ISOK    )   =  -1.0   ! ??
!!bc2xn_adv(IBC_ISOR2   ,IXADV_ISOR2   )   =  -1.0   ! SHort-lived
!!bc2xn_adv(IBC_HCOHCO  ,IXADV_HCOHCO  )   =  -1.0   ! Excluded
!!bc2xn_adv(IBC_CH3X    ,IXADV_CH3X    )   =  -1.0   ! ??
!!bc2xn_adv(IBC_NO3     ,IXADV_NO3     )   =  -1.0   ! SHort-lived
!!bc2xn_adv(IBC_N2O5    ,IXADV_N2O5    )   =   1.0   ! jej - 000927
!!bc2xn_adv(IBC_C3H7O2  ,IXADV_C3H7O2  )   =   1.0
!!bc2xn_adv(IBC_ACETONE ,IXADV_ACETON  )   =   1.0
!!bc2xn_adv(IBC_CH3COD  ,IXADV_CH3COD  )   =  -1.0   ! ??
!!bc2xn_adv(IBC_NOZ     ,IXADV_NOZ     )   =  -1.0   ! ???
!!bc2xn_adv(IBC_CH3O2   ,IXADV_CH3O2   )   =   Short-lived
!!bc2xn_adv(IBC_C2H5O2  ,IXADV_C2H5O2  )   =   Short-lived
!!bc2xn_adv(IBC_C4H9O2  ,IXADV_C4H9O2  )   =   Short-lived
!!bc2xn_adv(IBC_C6H13O2 ,IXADV_C6H13O2 )   =   Short-lived
!!bc2xn_adv(IBC_O3P     ,IXADV_O3P     )   =  -1.0   ! Short-lived
!!bc2xn_adv(IBC_O1D     ,IXADV_O1D     )   =  -1.0   ! Short-lived
!!bc2xn_adv(IBC_OH      ,IXADV_OH      )   =  -1.0   ! Short-lived
!!bc2xn_adv(IBC_O3NO    ,IXADV_O3NO    )   =  -1.0   ! Excluded
!!bc2xn_adv(IBC_DMS     ,IXADV_DMS     )   =  -1.0   ! Query ????
