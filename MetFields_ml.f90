module MetFields_ml
  use ModelConstants_ml,   only : KMAX_BND,KMAX_MID,NMET &
       ,IIFULLDOM, JJFULLDOM
  use Par_ml, only : MAXLIMAX,MAXLJMAX,GIMAX,GJMAX, me  &
       ,limax,ljmax,li0,li1,lj0,lj1 
  implicit none
  private


!----------------- basic met fields ----------------------------------!
!  Here we declare the meteorological fields used in the model        !
!
! Horizonal alignments....
! Placement of q(i,j), u(i,j), v(i,j) 
!
!    --------------- --------------- --------------- ---------------
!    |              |               |               |               | 
!    |              |               |               |               | 
!   u03   q13      u13    q23      u23     q33     u33    q43      u43 ... u(MAXLIMAX,3)
!    |              |               |               |               | 
!    |              |               |               |               | 
!    -----v12------- -----v22------- ------v32------ -----v42-------
!    |              |               |               |               | 
!    |              |               |               |               | 
!   u02   q12      u12    q22      u22     q32     u32    q42      u42
!    |              |               |               |               | 
!    |              |               |               |               | 
!    -----v11------- -----v21------- ------v31------ -----v41-------
!    |              |               |               |               | 
!    |              |               |               |               | 
!   u01   q11      u11    q21      u21     q31     u31    q41      u41
!    |              |               |               |               | 
!    |              |               |               |               | 
!    -----v10------- -----v20------- ------v30------ -----v40-------
!
  !---------------------------------------------------------------------!
  !
  !
  ! Vertical levels: z_mid,  z_bnd, sigma_mid, sigma_bnd
  !=============================================================================
  !*   "mid" and "bnd" are used as suffixes on z and sigma as shown in
  !*   the sketch below. "bnd" is the boundary between two layers and
  !*   "mid" the midddle of the layer. The numbering of layers starts
  !*   from 1 at the surface.
  !*
  !*
  !*
  !* ---------------------------
  !*
  !*
  !* - - - - - - - - - - - -     KMAX_MID -1
  !*
  !*
  !* --------------------------  KMAX_BDN-1       (z_bnd)   (sigma_bnd)
  !*
  !*
  !* - - - - - - - - -           KMAX_MID(old kmax2) = 20    (z_mid)   (sigma_mid)   (old z2)
  !*
  !* ------------------------    KMAX_BND = 21    (z_bnd)                 (old z1)
  !* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  surface \\\\\\\\\\\\\\\\
  !*
  ! RESULT FROM ONE TEST
  !Tested zm3d (=zm) which is sometimes used against
  !z_mid. Seems almost identical. Diff in exner functions maybe?
  !                   z_bnd        z_mid       zm3d 
  !     DEBUG_Z  2  1  16276.2359  15244.9023  15244.7131
  !     DEBUG_Z  2  2  14319.8682  13536.3354  13531.2981
  !     DEBUG_Z  2  3  12815.7707  12185.9090  12177.9740
  !     DEBUG_Z  2  4  11598.2202  11011.4342  11004.0286
  !     DEBUG_Z  2  5  10461.1309   9792.2354   9788.6193
  !     DEBUG_Z  2  6   9168.4701   8431.0622   8430.8546
  !     DEBUG_Z  2  7   7747.1534   7010.9878   7013.5860
  !     DEBUG_Z  2  8   6324.9278   5696.2857   5700.0226
  !     DEBUG_Z  2  9   5103.2253   4595.6436   4598.2893
  !     DEBUG_Z  2 10   4110.6562   3690.7686   3692.3391
  !     DEBUG_Z  2 11   3286.0770   2933.6493   2934.3594
  !     DEBUG_Z  2 12   2591.8972   2296.1632   2296.2730
  !     DEBUG_Z  2 13   2007.7924   1761.4943   1761.6936
  !     DEBUG_Z  2 14   1520.4434   1316.9528   1317.4565
  !     DEBUG_Z  2 15   1116.9477    950.9995    951.4678
  !     DEBUG_Z  2 16    787.3321    655.5900    655.8842
  !     DEBUG_Z  2 17    525.1633    424.5627    424.7443
  !     DEBUG_Z  2 18    324.8232    253.9725    254.0498
  !     DEBUG_Z  2 19    183.6019    137.4135    137.4284
  !     DEBUG_Z  2 20     91.3017     45.6234     45.6146


  !   Vertical level geopotential heights:

  real,public, save, &
       dimension(MAXLIMAX,MAXLJMAX,KMAX_BND) :: z_bnd ! height of full layers
  real,public, save, &
       dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: z_mid ! height of half layers

  !   Two sets of Met. fields are read in, and a linear interpolation is made
  !   between these two points in time. NMET  == 2 (two points in time)
  !   note u_xmj, v_xmi are not "real" m/s wind speeds
  !   - they are actually divided by the mapping factor in the perpendicular direction).
  !
  real,public, save, dimension(0:MAXLIMAX,MAXLJMAX,KMAX_MID,NMET) :: u_xmj 
  real,public, save, dimension(MAXLIMAX,0:MAXLJMAX,KMAX_MID,NMET) :: v_xmi


  real,public, save, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID,NMET) :: &
       th      &  ! Potential teperature  ( deg. k )
       ,q      &  ! Specific humidity
       ,roa    &  ! kg/m3
       ,cw        ! cloudwater
  real,public, save, dimension(MAXLIMAX,MAXLJMAX,KMAX_BND,NMET) :: &
        SigmaKz  &! vertical diffusivity in sigma coords
       ,sdot     &! vertical velocity, sigma coords, 1/s
       ,Kz_met    ! vertical diffusivity in sigma coordinates from meteorology


  ! since pr,cc3d,cc3dmax,cnvuf,cnvdf used only for 1 time layer - define without NMET
  real,public, save, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: &
        pr      & ! Precipitation
       ,cc3d    & ! 3-d cloud cover (cc3d),
       ,cc3dmax & ! and maximum for layers above a given layer
       ,lwc     & !liquid water content
  ! QUERY - should xksig be MID, not BND? Is it needed at all?
       ,Kz_m2s     ! estimated Kz, in intermediate sigma levels, m2/s

  real,public, save, dimension(MAXLIMAX,MAXLJMAX,KMAX_BND) :: &
        cnvuf   & ! convective_updraft_flux (kg/s/m2)
       ,cnvdf    ! convective_downdraft_flux (kg/s/m2)


 ! We don't need to calculate u,v for RiB, Kz for all layer in future maybe
 ! Still, for safety  we let this extent to K=1 for now

  real,public, save, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: &
        u_mid   & ! wind u-compnent, m/s (real, not projected)
       ,v_mid     ! wind v-compnent, m/s
  


! Surface fields, interpolated:
 real,public, save, dimension(MAXLIMAX,MAXLJMAX,NMET) :: &
        ps        &! Surface pressure Pa
       ,t2_nwp    & ! Temp 2 m   deg. K
       ,fh        & ! surf.flux.sens.heat W/m^2
       ,fl        & ! latent heat flux W/m^2
       ,tau       & ! surf. stress  N/m^2
  ! These fields only available for EMEP/PARLAM from 2002 on
       ,rh2m            & !  RH at 2m
       ,SoilWater       & !  Shallow  (Upper 7.2cm in PARLAM)
       ,SoilWater_deep  & !  Deep (Next 6x7cm in PARLAM), converted to relative value 
       ,sdepth          & !  Snowdepth, m
       ,ice_nwp             & ! QUERY why real?
       ,sst     &  ! SST Sea Surface Temprature- ONLY from 2002 in PARLAM
       ,ws_10m    ! wind speed 10m
 

 real,public, save, dimension(MAXLIMAX,MAXLJMAX) :: &
     u_ref             & ! wind speed m/s at 45m (real, not projected)
    ,rho_surf          & ! Surface density
    ,surface_precip    & ! Surface precip mm/hr
    ,Tpot2m            & ! Potential temp at 2m
    ,ustar_nwp         & ! friction velocity m/s ustar^2 = tau/roa
    ,invL_nwp          & ! friction velocity m/s ustar^2 = tau/roa
    ,pzpbl               ! stores H(ABL) for averaging and plotting purposes, m


!  temporary placement of solar radiation variations QUERY?
 
  real, public, dimension(MAXLIMAX, MAXLJMAX), save:: &
       zen          &  ! Zenith angle (degrees)
      ,coszen=0.0   &  ! cos of zenith angle
      ,Idiffuse     &  ! diffuse solar radiation (W/m^2)
      ,Idirect         ! total direct solar radiation (W/m^2)


 logical,public, save, dimension(MAXLIMAX,MAXLJMAX) :: &
       nwp_sea     ! Sea in NWP mode, determined in HIRLAM from roughness class

  real,public, save, dimension(MAXLIMAX,MAXLJMAX) :: &   !st-dust
       clay_frac  &  ! clay fraction (%) in the soil
      ,sand_frac     ! sand fraction (%) in the soil

  ! Different NWP outputs for soil water are possible. We can currently
  ! cope with two:
  character(len=10), public, save  :: SoilWaterSource  ! IFS or PARLAM

! Logical flags, used to determine if some met fields are present in the
! input or not:
  logical, public, save :: &
     foundustar     & ! Used for MM5-type, where u_xmj* but not tau
    ,foundsdot      & ! If not found: compute using divergence=0
    ,sdot_at_mid    & ! set false if sdot is defined
    ,foundSST       & ! false if no SeaSurfaceT in metdata
    ,foundSoilWater       & ! false if no SW-shallow
    ,foundSoilWater_deep  & ! false if no SW-deep
    ,foundsdepth    & ! false if no snow_flag depth in metdata
    ,foundice       & ! false if no ice_nwp coverage (%) in metdata
    ,foundnwp_sea   &  ! false if no rough file is found QUERY description?
  ! (when read) at level  boundaries and therefore do not need to be
  ! interpolated.
    ,foundKz_met    & ! false if no Kz from meteorology
    ,foundconv      & ! false if convection not found or not used
  ! Introduced for FUTURE NH3, but also sea-salt
    ,foundws10_met   & ! false if no u10 from meteorology
    ,foundu10_met   & ! false if no u10 from meteorology
    ,foundv10_met   & ! false if no v10 from meteorology
    ,foundprecip    & ! false if no precipitationfrom meteorology
    ,foundcloudwater  !false if no cloudwater found

end module MetFields_ml
