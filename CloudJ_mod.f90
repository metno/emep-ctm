  !-----------------------------------------------------------------------!
  !
  !     This is the main module that calculates J-values in the EMEP model,
  !     using the CloudJ modules stored in ModsCloudJ_mod.f90.
  !    
  !-----------------------------------------------------------------------!

MODULE CloudJ_mod
  
    use FJX_CMN_MOD  
    use FJX_SUB_MOD
    use FJX_INIT_MOD 
    use CLD_SUB_MOD,           only: CLOUD_JX

    use DefPhotolysis_mod      
    use GridValues_mod,        only: glat,glon,A_bnd,B_bnd,A_mid,B_mid,dA,dB
    use LandDefs_mod,          only: LandType, LandDefs
    use Landuse_mod,           only: LandCover
    use MetFields_mod,         only: ps,foundcloudwater,q,th,cc3d,  &
                                     roa, z_bnd, cw_met, &
                                     foundcloudicewater, ciw_met
    use GasParticleCoeffs_mod, only: DDspec
    use Config_module,         only: KMAX_BND,KMAX_MID,KCHEMTOP,METSTEP,USES, &
                                     NPROC, IOU_INST, num_lev3d,lev3d, cloudjx_strat, &
                                     MasterProc
    use Par_mod,               only: me,LIMAX,LJMAX
    use TimeDate_mod,          only: daynumber,current_date,date
    use ZchemData_mod,         only: rcphot, rcphotslice, rh
    use Functions_mod,         only: Tpot_2_T
    use TimeDate_ExtraUtil_mod,only: date2string
    use SmallUtils_mod,        only: find_index
    use DerivedFields_mod,     only: d_3d, f_3d

    use Chemfields_mod,        only: xn_adv, xn_bgn, xn_shl, NSPEC_BGN, Dobson
    use ChemDims_mod,          only: NSPEC_ADV, NPHOTOLRATES
    use ChemSpecs_mod
    use ChemRates_mod,         only: photol_used, setPhotolUsed
    use PhysicalConstants_mod, only: AVOG, ATWAIR
    use AeroFunctions_mod,     only: LogNormFracBelow, GerberWetRad, GerberWetSig, pmH2O_gerberSig, DpgN2DpgV
    use CheckStop_mod,         only: StopAll
                                  

    IMPLICIT NONE
    PUBLIC

    real, parameter ::  PI180 = 3.141592653589793d0/180.d0 
    real, parameter ::  molec_cm3 = 6.022140857e17 ! convert mol/m3 to molec/cm3
    real, parameter ::  dobson_km = 2.24115e06     ! convert mol/m3 to Dobson Unit/km
    real  GMTAU, ALBEDO, XGRD, YGRD, PSURF, SCALEH
    real  CF, PMID, ZDEL, ICWC, F1
    integer I,J,K,L,N
    integer NRAN, NJXX

!     --------------------params sent to CLOUD_JX-------------------------
    real                       :: U0,SZA,REFLB,SOLF,CLDCOR,FG0
    logical, save              :: LPRTJ=.false.
    integer                    :: CLDFLAG,NRANDO,IRAN,LNRG,ICNT
    integer                    :: NICA,JCOUNT
    character*6, dimension(JVN_)  ::  TITLJXX

    integer, save :: photo_out_ix_no2_fj    = -1
    integer, save :: photo_out_ix_o3a_fj    = -1
    integer, save :: photo_out_ix_o3b_fj    = -1
    integer, save :: photo_out_ix_h2o2_fj   = -1
    integer, save :: photo_out_ix_hno3_fj   = -1
    integer, save :: photo_out_ix_ach2o_fj  = -1
    integer, save :: photo_out_ix_bch2o_fj  = -1
    integer, save :: photo_out_ix_hono_fj   = -1
    integer, save :: photo_out_ix_ho2no2_fj = -1
    integer, save :: photo_out_ix_no3_fj    = -1
    integer, save :: photo_out_ix_ch3o2h_fj = -1 
    integer, save :: photo_out_ix_MEK_fj    = -1
    integer, save :: photo_out_ix_N2O5_fj   = -1
    integer, save :: photo_out_ix_GLYOX_fj  = -1 ! HCOHCO 
    integer, save :: photo_out_ix_CH3CHO_fj = -1 
    integer, save :: photo_out_ix_ACETON_fj = -1
    integer, save :: photo_out_ix_MGLYOX_fj = -1 ! RCOCHO
    integer, save :: photo_out_ix_BIACET_fj = -1 ! IDCH3COY
    integer, save :: photo_out_ix_PAN_fj    = -1 
    integer, save :: photo_out_ix_GLYOXA_fj = -1
    integer, save :: photo_out_ix_GLYOXB_fj = -1 
    integer, save :: photo_out_ix_GLYOXC_fj = -1 

    integer, save :: sulph_i, dustfroad_i, dustfwb_i, dustfsah_i
    integer, save :: dustcroad_i, dustcwb_i, dustcsah_i
    integer, save :: seasaltf_i, seasaltc_i
    integer, save :: ffirebc_i, ffireom_i, ffireppm25_i, ffirec_i 

    integer, save :: DDdep_ssf, DDdep_ssc, DDdep_dustf, DDdep_dustc
    
!---U0 = cos (SZA), SZA = solar zenith angle. Activate GIT?
!---REFLB = Lambertian reflectivity at the Lower Boundary
!---SOLF = solar flux factor for sun-earth distance
!---FG0 = scale for asymmetry factor to get equivalent isotropic (CLDFLAG=3 only)
!---LPRTJ = .true. = turn on internal print in both CLOUD_JX & PHOTO_JX
!--- P = edge press (hPa), Z = edge alt (m), T = layer temp (K)
!--- D = layer dens-path (# molec /cm2), O = layer O3 path (# O3 /cm2)
!--- R = layer rel.hum.(fraction)
!---LWP/IWP = Liquid/Ice water path (g/m2)
!---REFFL/REFFI = R-effective(microns) in liquid/ice cloud
!---CLF = cloud fraction (0.0 to 1.0)
!---CLDIW = integer denoting cloud in layer: 1=water, 2=ice, 3=both
!---AERSP = aerosol path (g/m2) & NDXAER = aerosol index type
!---  aerosols are dimensioned with up to AN_ different types in an ICA layer
!---L1_ = parameter, dim of profile variables, L_+1 for top (non CTM) layer
!---AN_ = parameter, dim of number of aerosols being passed (25 standard)
!---VALJXX = J-values from CLOUD_JX & PHOTO_JX
!---JVN_ = dim of max number of J-s reported out (in the order of fast-JX, not CTM)
!---CLDFLAG = integer index for type of cloud overlap
!---CLOUD_JX:   different cloud schemes
!---CLOUD_JX:   different cloud schemes (4:8 require max-ran overlap algorithm)
!       CLDFLAG = 1  :  Clear sky J's
!       CLDFLAG = 2  :  Averaged cloud cover
!       CLDFLAG = 3  :  cloud-fract**3/2, then average cloud cover
!       CLDFLAG = 4  :  Average direct solar beam over all ICAs, invert to get clouds
!       CLDFLAG = 5  :  Random select NRANDO ICA's from all(Independent Column Atmos.)
!       CLDFLAG = 6  :  Use all (up to 4) quadrature cloud cover QCAs (mid-pts of bin)
!       CLDFLAG = 7  :  Use all (up to 4) QCAs (average clouds within each Q-bin)
!       CLDFLAG = 8  :  Calculate J's for ALL ICAs (up to 20,000 per cell!)
!---NRANDO = number of random ICAs to do J's for (CLDFLAG=4)
!---IRAN = starting index for random number selection
!---LNRG = flag for setting max-ran overlap groups:
!---     =0   break max overlap groups at cloud fraction = 0
!---     =3   else fixed 3 layers (1:9, 9:last LWcloud, LWclud+1:LTOP)
!---     else(=6) fixed correlated length max-overlap layers
!---NICA = total number of ICAs

!-----these are the key results coming out of fast-JX core
!-----they need to be dimensioned here and passed to fast-JX to fill.

CONTAINS

SUBROUTINE setup_phot_cloudj(i_emep,j_emep,errcode,mode)

    !calculate rcphot for one ICA in the emep model
    IMPLICIT NONE

    integer, intent(in) :: i_emep, j_emep
    integer, intent(inout) :: errcode
    integer,  intent(in) :: mode
    integer :: ilu, lu, nr_local
    real :: Pres_mid, temperature, swp

    integer IDAY, JLAT, ILON

    real, allocatable, save :: CWIC(:), CWLC(:)
    real, allocatable, save :: DeltZ(:),CLWP(:),CIWP(:)
    real, allocatable, save :: AER(:,:)
    integer, allocatable, save :: NAA(:,:)
    real, allocatable, save :: DUST_F(:),DUST_C(:),SULF(:),SEAS_F(:),SEAS_C(:),Biomas_B(:),Forest_F(:)

    logical, save :: first_call=.true.

    ! ICA input to CloudJ. NB: L_in sets the number of levels read in by CLOUD_JX
    real, allocatable, save :: PPP(:), ZZZ(:), TTT(:), DDD(:), RRR(:)
    real, allocatable, save :: REFFI(:), CLF(:), OOO(:), LWP(:), REFFL(:), IWP(:)
    real, allocatable, save :: AERSP(:,:)

    integer, allocatable, save :: CLDIW(:)
    integer, allocatable, save :: NDXAER(:,:)
    real,  allocatable, save :: VALJXX(:,:)

    real :: BB_BCratio_real, FF_BCratio_real
    integer :: BB_BC_index, FF_BC_index

    ! declarations used in reading/handling satellite stratospheric ozone data
    integer :: la, lo, w, z
    integer, save :: dim_lon, dim_lat, dim_alt, OZ_TOP, oz_month=-999, IO_OZONE=78, NAERMAX=25
    integer, save :: year, month, naerosol=15 ! max number of included EMEP aerosol
    real, allocatable, save :: lon_ozone(:), lat_ozone(:), hpa_ozone(:)
    real, allocatable, save :: temp_ozone(:,:,:), pres_ozone(:,:,:), ozon_ozone(:,:,:)
    real, allocatable, save :: temp_obs(:,:,:), pres_obs(:,:,:), ozon_obs(:,:,:)
    character(len=150), save  :: fname_ozone
    real,  save :: dlat, dlon, emep_top
    integer, save :: satellite_altind

    ! mass fraction mapping between EMEP fine/coarse mode aerosol and CloudJ UMich bins
    real, save  :: SS_bin1, SS_bin2, SS_bin3 ,SS_bin4 
    real, save  :: DU_bin1, DU_bin2, DU_bin3 ,DU_bin4
    real, save  :: dry_d, wet_d_ssf, sig_wet_ssf, wet_d_ssc, sig_wet_ssc, pmH2O_wet
    real, save  :: wet_mmd_ssf, wet_mmd_ssc, dust_f_mmd, dust_c_mmd
    
    !---fast-JX:  INIT_JX is called only once to read in & store all fast-JX data: 
    !             also sets up random sequence for cloud-JX. 
    ! JVN_ = 101 (max number of J values read in from CloudJ)
    ! NJXX = number of derived Jvalues, set in RD_XXX 

    ! set time-variables; nr_local = 1 uses current meteorology timestep 
    nr_local=1
    if(mode>0)nr_local=mode
    iday = daynumber  
    year = current_date%year
    month = current_date%month
    gmtau = current_date%hour + current_date%seconds/3600.0
    if(nr_local==2)gmtau = gmtau + METSTEP ! use next met tstep

    ! read monthly lat-lon-altitude ozone satellite observations
    if(oz_month/=month) then 
          ! ##########################################################################################
          ! Ozone satellite observation files are preprocessed using the script Sat_Ozone_datsaver.py.
          ! Data is defined for altitude levels strictly above 100 hPa, so as not to                  
          ! conflict with the model top of the (current) EMEP model.                                  
          !                                                                                                            
          ! Data is based on measurements from six satellite missions between yrs 2005-2021:   
          !     https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-ozone-v1?tab=form
          !
          ! ##########################################################################################

          ! define path based on model year & month and read in monthly lon-lat-alt satellite obs. data        
          if(year < 2005 .or. year > 2021 .or. USES%CLIMSTRATO3) then 
            fname_ozone = trim(cloudjx_strat)//date2string("/clim_MM.dat",current_date)
          else
            fname_ozone = trim(cloudjx_strat)//date2string("/YYYYMM.dat",current_date)
          endif
          if(me==0)write(*,*) 'Opening satellite O3/T obs. file: ', fname_ozone 
          
          open(IO_OZONE, file=fname_ozone,form='unformatted')
          read(IO_OZONE) dim_lon, dim_lat, dim_alt

          if(first_call) then 
                ! i/o arrays necessary for strat ozone files
                allocate(lon_ozone(dim_lon))
                allocate(lat_ozone(dim_lat))
                allocate(hpa_ozone(dim_alt))

                allocate(temp_ozone(dim_lon,dim_lat,dim_alt))
                allocate(ozon_ozone(dim_lon,dim_lat,dim_alt))

                ! arrays to hold the satellite data; pres at mid-lvls, temp and O3 at half-lvls (starting below 1st mid p-lvl)
                allocate(temp_obs(LIMAX,LJMAX,dim_alt))
                allocate(pres_obs(LIMAX,LJMAX,dim_alt))
                allocate(ozon_obs(LIMAX,LJMAX,dim_alt))
          end if

          read(IO_OZONE) lon_ozone ! -180 to 180 degrees, as in EMEP
          read(IO_OZONE) lat_ozone ! -90 to 90 degrees, as in EMEP
          read(IO_OZONE) hpa_ozone ! interface levels in hPa 

          dlon = lon_ozone(2) - lon_ozone(1) ! grid spacing in degrees
          dlat = lat_ozone(2) - lat_ozone(1) 

          read(IO_OZONE) temp_ozone ! temperature (K) mid-lvls
          read(IO_OZONE) ozon_ozone ! ozone (molec/cm2) mid-lvls

          close(IO_OZONE)     

          do w=1,LIMAX ! map lat(LJMAX) and lon(LIMAX) indices and populate the obs arrays
                do z=1,LJMAX 
                      lo = max(1,int( 1 / dlon * glon(w,z) + 179.999 / dlon + 1))
                      la = max(1,int( 1 / dlat * glat(w,z) + 89.999 / dlat + 1))

                      lo = min(dim_lon, lo) ! avoid any overshoot (e.g. when lon = 180.25 deg)
                      la = min(dim_lat, la)

                      temp_obs(w,z,:) = temp_ozone(lo,la,:)
                      ozon_obs(w,z,:) = ozon_ozone(lo,la,:)
                end do
          end do

          if (first_call) then
                ! calculate startpoint satellite lvls above emep model top
                emep_top = (A_bnd(1) + B_bnd(1)*1e5) / 1e2 ! hPa based on mean surf_P; might be lower in some places
                satellite_altind = 1

                ! for information; satellite data starts at 200 hPa 
                do while (emep_top < hpa_ozone(satellite_altind)) 
                  satellite_altind = satellite_altind + 1
                end do

                ! pressure in Pa here
                if(MasterProc .and. USES%EtaCOORDINATES .and. A_bnd(1) + B_bnd(1)*1e5 > 2.e4) write(*,*) & 
                  'Warning: CloudJ satellite O3/T constant between 200 hPa down to model top.'

                ! total number of levels with O3/T data stacked on top EMEP levels
                OZ_TOP = KMAX_MID + dim_alt - satellite_altind

                ! allocate i/o arrays necessary for cloudj
                allocate(PPP(OZ_TOP+1), ZZZ(OZ_TOP+1)) 
                allocate(TTT(OZ_TOP),   DDD(OZ_TOP),   RRR(OZ_TOP),   CLF(OZ_TOP))
                allocate(OOO(OZ_TOP),   LWP(OZ_TOP),   IWP(OZ_TOP)               )
                allocate(REFFI(OZ_TOP), CLDIW(OZ_TOP), REFFL(OZ_TOP)             )
                allocate(NDXAER(OZ_TOP,naerosol),      AERSP(OZ_TOP,naerosol)    )
                allocate(VALJXX(OZ_TOP-1,JVN_)                                   )

                allocate(CWIC(OZ_TOP),         CWLC(OZ_TOP),   DeltZ(OZ_TOP), CLWP(OZ_TOP))
                allocate(DUST_F(OZ_TOP),       DUST_C(OZ_TOP), CIWP(OZ_TOP)               )
                allocate(SULF(OZ_TOP),         SEAS_F(OZ_TOP), SEAS_C(OZ_TOP)             )
                allocate(Biomas_B(OZ_TOP),     Forest_F(OZ_TOP)                           )
                allocate(NAA(OZ_TOP,naerosol), AER(OZ_TOP,naerosol)                       )
          endif

    end if ! oz_month /= month
    
    oz_month = month

    ! initialize requested (in config file) J-value output arrays & ModsCloudJ config
    if(first_call)then 
          ! parameters from FJX_CMN_MOD used to initialize cloud-j arrays
          LPAR   = OZ_TOP              ! this can be set by CTM code. 
          LWEPAR = KMAX_MID - KCHEMTOP ! LWEPAR = number of lvls with cloud calcs
          L_     = LPAR
          L1_    = L_ + 1              ! L_ = number of CTM layers
          L2_    = 2*L_ + 2            ! no. levels in the Fast-JX grid that
                                       ! includes both layer edges and layer mid-points
          JVL_   = KMAX_MID - KCHEMTOP ! vertical(levels) dim for J-values sent to CTM
          
          call INIT_FJX(TITLJXX,JVN_,NJXX) 

          if(allocated(f_3d))then ! if output fields requested in config file, set index
            photo_out_ix_no2_fj    = find_index("D3_J(NO2)",    f_3d(:)%subclass)
            photo_out_ix_o3a_fj    = find_index("D3_J(O3a)",    f_3d(:)%subclass)
            photo_out_ix_o3b_fj    = find_index("D3_J(O3b)",    f_3d(:)%subclass)
            photo_out_ix_h2o2_fj   = find_index("D3_J(H2O2)",   f_3d(:)%subclass)
            photo_out_ix_hno3_fj   = find_index("D3_J(HNO3)",   f_3d(:)%subclass)
            photo_out_ix_ach2o_fj  = find_index("D3_J(ACH2O)",  f_3d(:)%subclass)
            photo_out_ix_bch2o_fj  = find_index("D3_J(BCH2O)",  f_3d(:)%subclass)
            photo_out_ix_hono_fj   = find_index("D3_J(HONO)",   f_3d(:)%subclass)
            photo_out_ix_ho2no2_fj = find_index("D3_J(HO2NO2)", f_3d(:)%subclass)
            photo_out_ix_no3_fj    = find_index("D3_J(NO3)",    f_3d(:)%subclass)
            photo_out_ix_ch3o2h_fj = find_index("D3_J(CH3O2H)", f_3d(:)%subclass)
            photo_out_ix_MEK_fj    = find_index("D3_J(MEK)",    f_3d(:)%subclass)
            photo_out_ix_N2O5_fj   = find_index("D3_J(N2O5)",   f_3d(:)%subclass)
            photo_out_ix_GLYOX_fj  = find_index("D3_J(GLYOX)",  f_3d(:)%subclass)
            photo_out_ix_CH3CHO_fj = find_index("D3_J(CH3CHO)", f_3d(:)%subclass)
            photo_out_ix_ACETON_fj = find_index("D3_J(ACETON)", f_3d(:)%subclass)
            photo_out_ix_MGLYOX_fj = find_index("D3_J(MGLYOX)", f_3d(:)%subclass)
            photo_out_ix_BIACET_fj = find_index("D3_J(BIACET)", f_3d(:)%subclass)
            photo_out_ix_PAN_fj    = find_index("D3_J(PAN)",    f_3d(:)%subclass)
            photo_out_ix_GLYOXA_fj = find_index("D3_J(GLYOXA)", f_3d(:)%subclass)
            photo_out_ix_GLYOXB_fj = find_index("D3_J(GLYOXB)", f_3d(:)%subclass)
            photo_out_ix_GLYOXC_fj = find_index("D3_J(GLYOXC)", f_3d(:)%subclass)
            if(MasterProc)write(*,*) 'Outputting CloudJ J-values specified in config.'
          endif

          if (MasterProc) write(*,*) 'Cloud water (1) and ice (2) for photolysis rate calculations: ' &
                     ,foundcloudwater, foundcloudicewater

          if(MasterProc .and. .not. foundcloudicewater) write(*,*) 'WARNING: Running CloudJ without cloud ice water content: ' &
                                                , 'Might overestimate surface photolysis rates.'
    endif

    ! find available aerosol indices from chemistry
    if(first_call)then

      sulph_i      = find_index( "SO4           ", species_adv(:)%name )
      dustfroad_i  = find_index( "Dust_road_f   ", species_adv(:)%name )
      dustfwb_i    = find_index( "Dust_wb_f     ", species_adv(:)%name )
      dustfsah_i   = find_index( "Dust_sah_f    ", species_adv(:)%name )
      dustcroad_i  = find_index( "Dust_road_c   ", species_adv(:)%name )
      dustcwb_i    = find_index( "Dust_wb_c     ", species_adv(:)%name )
      dustcsah_i   = find_index( "Dust_sah_c    ", species_adv(:)%name )
      seasaltf_i   = find_index( "SeaSalt_f     ", species_adv(:)%name )
      seasaltc_i   = find_index( "SeaSalt_c     ", species_adv(:)%name )
      ffirebc_i    = find_index( "ffire_BC      ", species_adv(:)%name )
      ffireom_i    = find_index( "ffire_OM      ", species_adv(:)%name )
      ffireppm25_i = find_index( "ffire_remPPM25", species_adv(:)%name )
      ffirec_i     = find_index( "ffire_c       ", species_adv(:)%name )

      ! aerosol properties 
      DDdep_ssf   = find_index(  'SSf',DDspec(:)%name, any_case=.true. )
      DDdep_ssc   = find_index(  'SSc',DDspec(:)%name, any_case=.true. )
      DDdep_dustf = find_index(  'DUf',DDspec(:)%name, any_case=.true. )
      DDdep_dustc = find_index(  'DUc',DDspec(:)%name, any_case=.true. )

    endif

    !inputs for SOLAR_JX 
    YGRD = glat(i_emep,j_emep)*PI180
    XGRD = glon(i_emep,j_emep)*PI180
    PSURF = ps(i_emep,j_emep,nr_local)/100.0 !Pa-> hPa

    ILON = int(XGRD/PI180)
    JLAT = int(YGRD/PI180)
    IRAN = 13+ILON+3*JLAT+5*(YEAR-1900)+7*IDAY + 11*nint(GMTAU) ! random number based on year/day/hour

    ! initialize cloud fraction, water path and aerosol arrays to zero
    CLF  = 0.
    CLWP = 0.
    CIWP = 0.
    CWIC = 0.
    CWLC = 0.
    RRR  = 0.

    SULF = 0.
    DUST_F = 0.
    DUST_C = 0.
    SEAS_F = 0.
    SEAS_C = 0.
    Forest_F = 0.
    Biomas_B = 0.

    AER = 0.
    NAA = 0

    !use fastj vertical direction, i.e. L largest at top, 1 at surface 
    L=0

    ! EMEP mid-levels 
    do k=kmax_mid,1,-1
          L = L+1
          Pres_mid = A_mid(k) + B_mid(k)*ps(i_emep,j_emep,nr_local)
    
          !potential -> absolute temperature
          temperature = th(i_emep,j_emep,k,nr_local)* Tpot_2_T( Pres_mid )
          TTT(L) = temperature

          if (k >= KCHEMTOP) RRR(L) = rh(k) ! use RH calculated in setup_1d; only defined up to chemtop

          !cloud cover stored as fraction from 0 to 1 (line 652 Met_mod.f90)
          DeltZ(L) = z_bnd(i_emep,j_emep,k)-z_bnd(i_emep,j_emep,k+1)
          CF = cc3d(i_emep,j_emep,k,1)

          !if clouds are present, populate cloud liquid and ice water content arrays
          if(CF > 0.00001d0) then
                CLF(L) = CF
                if(foundcloudwater) then
                   CWLC(L) = cw_met(i_emep,j_emep,k,1) / CLF(L)  ! kg/kg per cloud area 
                   CLWP(L) = CWLC(L) * 1000. * roa(i_emep,j_emep,k,1) * DeltZ(L)  ! water path [g/m2]
                endif  
                if(foundcloudicewater) then
                   CWIC(L) = ciw_met(i_emep,j_emep,k,1) / CLF(L) ! kg/kg per cloud area 
                   CIWP(L) = CWIC(L) * 1000. * roa(i_emep,j_emep,k,1) * DeltZ(L)  ! water path [g/m2]
                endif
          end if

          ! negative aerosol indices give UMich aerosol definitions (include RH effects). FF and BB as a function of 
          ! black carbon (BC) mass content (%). FF = Fossil Fuel, BB = Biomass Burning.
          ! 
          ! 1 SULF    2 SS-1    3 SS-2    4 SS-3    5 SS-4    6 DD-1    7 DD-2
          ! 8 DD-3    9 DD-4    10 FF00   11 FF02   12 FF04   13 FF06   14 FF08
          ! 15 FF10   16 FF12   17 FF14   18 BB00   19 BB02   20 BB04   21 BB06
          ! 22 BB08   23 BB10   24 BB12   25 BB14   26 BB16   27 BB18   28 BB20
          ! 29 BB22   30 BB24   31 BB26   32 BB28   33 BB30

          IF(USES%CLOUDJAEROSOL) then
            ! emissions are in molecules / cm3 / s and are used to populate xn_2d,
            ! which is then converted to ppb before being put into xn_adv (i.e., a ratio).
            !
            ! atmos. densities from kg/m3 to g/m3 to moles/m3 to moles/m2 to grams/m2 of aerosol
            if(sulph_i>0) SULF(L)   = xn_adv(sulph_i,i_emep,j_emep,k)            &
                                    * 1000 * roa(i_emep,j_emep,k,1) / ATWAIR     &
                                    * DeltZ(L) * species_adv(sulph_i)%molwt

            ! ----fine mode dust 
            if(dustfroad_i>0) DUST_F(L) = xn_adv(dustfroad_i,i_emep,j_emep,k)              &
                                        * 1000 * roa(i_emep,j_emep,k,1) / ATWAIR           &
                                        * DeltZ(L) * species_adv(dustfroad_i)%molwt

            if(dustfwb_i>0) DUST_F(L)   = DUST_F(L) + xn_adv(dustfwb_i,i_emep,j_emep,k)    & 
                                        * 1000 * roa(i_emep,j_emep,k,1) / ATWAIR           &
                                        * DeltZ(L) * species_adv(dustfwb_i)%molwt

            if(dustfsah_i>0) DUST_F(L)  = DUST_F(L) + xn_adv(dustfsah_i,i_emep,j_emep,k)   & 
                                        * 1000 * roa(i_emep,j_emep,k,1) / ATWAIR           &
                                        * DeltZ(L) * species_adv(dustfsah_i)%molwt

            ! ----coarse mode dust
            if(dustcroad_i>0) DUST_C(L) = xn_adv(dustcroad_i,i_emep,j_emep,k)              & 
                                        * 1000 * roa(i_emep,j_emep,k,1) / ATWAIR           &
                                        * DeltZ(L) * species_adv(dustcroad_i)%molwt

            if(dustcwb_i>0)   DUST_C(L) = DUST_C(L) + xn_adv(dustcwb_i,i_emep,j_emep,k)    & 
                                        * 1000 * roa(i_emep,j_emep,k,1) / ATWAIR           &
                                        * DeltZ(L) * species_adv(dustcwb_i)%molwt

            if(dustcsah_i>0)  DUST_C(L) = DUST_C(L) + xn_adv(dustcsah_i,i_emep,j_emep,k)   & 
                                        * 1000 * roa(i_emep,j_emep,k,1) / ATWAIR           &
                                        * DeltZ(L) * species_adv(dustcsah_i)%molwt

            ! ----fine sea salt; hygroscopic growth mass gain using Gerber functions
            if (seasaltf_i>0) then 

              SEAS_F(L)  = xn_adv(seasaltf_i,i_emep,j_emep,k)               & 
                         * 1000 * roa(i_emep,j_emep,k,1) / ATWAIR           &
                         * species_adv(seasaltf_i)%molwt   ! here still in g/m3

              ! calculate H2O mass using Gerber functions for geom. mean radius growth and geom. std growth
              dry_d       = DDspec(DDdep_ssf)%DpgN ! geometric mean diameter in meters
              wet_d_ssf   = 2 * GerberWetRad( dry_d * 0.5, RRR(L), DDspec(DDdep_ssf)%Gb )
              sig_wet_ssf = GerberWetSig( dry_d * 0.5, DDspec(DDdep_ssf)%sigma, RRR(L), DDspec(DDdep_ssf)%Gb )

              wet_mmd_ssf = DpgN2DpgV( wet_d_ssf, sig_wet_ssf ) ! mass median diam. for use in mass-frac calculations

              pmH2O_wet = pmH2O_gerberSig( SEAS_F(L) * 1e6,                          &  ! micrograms of seasalt per m3 as input
                          rho_kgm3=DDspec(DDdep_ssf)%rho_p,Dp=dry_d,Dpw=wet_d_ssf,   &
                          sigma_d=DDspec(DDdep_ssf)%sigma,sigma_w=sig_wet_ssf ) * 1e-6  ! output in micrograms per m3, convert to g/m3

              SEAS_F(L) = ( SEAS_F(L) + pmH2O_wet ) * DeltZ(L) ! g/m2 of wet fine sea salt

            endif ! fine sea salt
            
            ! ----coarse sea salt; include hygroscopic growth effect using Gerber functions
            if (seasaltc_i>0) then

              SEAS_C(L)  = xn_adv(seasaltc_i,i_emep,j_emep,k)               &
                         * 1000 * roa(i_emep,j_emep,k,1) / ATWAIR           &
                         * species_adv(seasaltc_i)%molwt   ! here still in g/m3

              ! calculate H2O mass using Gerber functions for geom. mean radius growth and geom. std growth
              dry_d       = DDspec(DDdep_ssc)%DpgN ! geometric mean diameter in meters
              wet_d_ssc   = 2 * GerberWetRad( dry_d * 0.5, RRR(L), DDspec(DDdep_ssc)%Gb )
              sig_wet_ssc = GerberWetSig( dry_d * 0.5, DDspec(DDdep_ssc)%sigma, RRR(L), DDspec(DDdep_ssc)%Gb )

              wet_mmd_ssc = DpgN2DpgV( wet_d_ssc, sig_wet_ssc ) ! mass median diam. for use in mass-frac calculations

              pmH2O_wet = pmH2O_gerberSig( SEAS_C(L) * 1e6,                          &  ! micrograms of seasalt per m3 as input
                          rho_kgm3=DDspec(DDdep_ssc)%rho_p,Dp=dry_d,Dpw=wet_d_ssc,   &
                          sigma_d=DDspec(DDdep_ssc)%sigma,sigma_w=sig_wet_ssc ) * 1e-6  ! output in micrograms per m3, convert to g/m3

              SEAS_C(L) = ( SEAS_C(L) + pmH2O_wet ) * DeltZ(L) ! g/m2 of wet coarse sea salt

            endif ! coarse sea salt
     
            ! ----forest fire
            if(ffirebc_i>0) Forest_F(L)   = xn_adv(ffirebc_i,i_emep,j_emep,k)                 & 
                                          * 1000 * roa(i_emep,j_emep,k,1) / ATWAIR            &
                                          * DeltZ(L) * species_adv(ffirebc_i)%molwt

            ! CMXbb:# The tricky bit: we read in BC, OC and PM25 from FINN, but want BC, OM and
            ! CMXbb:# rempPM25 for EmChem19a , CRI etc. We solve this by using FINN's
            ! CMXbb:# OC to estimate OM (factor 1.7), and subtracting both BC and OM from
            ! CMXbb:# PM25 to get remPM25.   (ForestFire_mod will prevent zeros).
            !
            ! hence remppm2.5 = pm2.5 - 1.7*OC - 1.0*BC
            ! such that for cloudj, forest fire pm2.5 = remppm2.5 + 1.7*OC (=OM) + 1.0*BC
            if(ffireom_i>0) Forest_F(L)   = Forest_F(L) + xn_adv(ffireom_i,i_emep,j_emep,k)   & 
                                          * 1000 * roa(i_emep,j_emep,k,1) / ATWAIR            &
                                          * DeltZ(L) * species_adv(ffireom_i)%molwt

            if(ffireppm25_i>0) Forest_F(L)= Forest_F(L) + xn_adv(ffireppm25_i,i_emep,j_emep,k)& 
                                          * 1000 * roa(i_emep,j_emep,k,1) / ATWAIR            &
                                          * DeltZ(L) * species_adv(ffireppm25_i)%molwt

            ! GFAS forest fire can also include a coarse particle mode
            if(ffirec_i>0) Forest_F(L)   = Forest_F(L) + xn_adv(ffirec_i,i_emep,j_emep,k)     & 
                                          * 1000 * roa(i_emep,j_emep,k,1) / ATWAIR            &
                                          * DeltZ(L) * species_adv(ffirec_i)%molwt
           
            ! naerosol variable must be set with the number of used aeorsol species
            AER(L,1)=SULF(L)   ! aerosol [g/m2]
            NAA(L,1)=-1        ! aerosol index

            ! UMich sea salt bin RADIUS range 0.05-0.63, 0.63-1.26, 1.26-2.50, 2.50-10 um
            ! EMEP fine mode and coarse mode from wet radius calculations; bin1 covers basically all of fine mode
            if (seasaltf_i>0) SS_bin1 = (LogNormFracBelow(wet_mmd_ssf,sig_wet_ssf,1.26e-6) &
                                       - LogNormFracBelow(wet_mmd_ssf,sig_wet_ssf,0.10e-6) )
        
            if (seasaltc_i>0) then 

              SS_bin2 = (LogNormFracBelow(wet_mmd_ssc,sig_wet_ssc,2.52e-6) &
                       - LogNormFracBelow(wet_mmd_ssc,sig_wet_ssc,1.26e-6) ) 

              SS_bin3 = (LogNormFracBelow(wet_mmd_ssc,sig_wet_ssc,5.00e-6) & 
                       - LogNormFracBelow(wet_mmd_ssc,sig_wet_ssc,2.52e-6) )

              SS_bin4 = (LogNormFracBelow(wet_mmd_ssc,sig_wet_ssc,20.00e-6)&
                       - LogNormFracBelow(wet_mmd_ssc,sig_wet_ssc,5.00e-6) )
                       
            end if

            AER(L,2)=SEAS_F(L) * SS_bin1 ! aerosol [g/m2] 
            NAA(L,2)=-2                  ! aerosol index

            AER(L,3)=SEAS_C(L) * SS_bin2 ! aerosol [g/m2] 
            NAA(L,3)=-3                  ! aerosol index

            AER(L,4)=SEAS_C(L) * SS_bin3 ! aerosol [g/m2] 
            NAA(L,4)=-4                  ! aerosol index

            AER(L,5)=SEAS_C(L) * SS_bin4 ! aerosol [g/m2]
            NAA(L,5)=-5                  ! aerosol index   

            ! UMich dust bin RADIUS range  
            ! 0.05-0.63, 0.63-1.26, 1.26-2.50, 2.50-10 um
            ! EMEP fine mode dust MMDiameter = 0.33 um, sigma = 1.8 
            ! EMEP coarse mode dust MMDiameter = 5.00 um, sigma = 2.2
            if(first_call .and. L==1) then ! dust bins are time-invariant because no hygroscopic growth

              dust_f_mmd = DpgN2DpgV( DDspec(DDdep_dustf)%DpgN, DDspec(DDdep_dustf)%sigma )
              dust_c_mmd = DpgN2DpgV( DDspec(DDdep_dustc)%DpgN, DDspec(DDdep_dustc)%sigma )

              DU_bin1 = (LogNormFracBelow(dust_f_mmd,DDspec(DDdep_dustf)%sigma,1.26e-6) &
                       - LogNormFracBelow(dust_f_mmd,DDspec(DDdep_dustf)%sigma,0.10e-6) )
          
              DU_bin2 = (LogNormFracBelow(dust_c_mmd,DDspec(DDdep_dustc)%sigma,2.52e-6) &
                       - LogNormFracBelow(dust_c_mmd,DDspec(DDdep_dustc)%sigma,1.26e-6) )

              DU_bin3 = (LogNormFracBelow(dust_c_mmd,DDspec(DDdep_dustc)%sigma,5.00e-6) &
                       - LogNormFracBelow(dust_c_mmd,DDspec(DDdep_dustc)%sigma,2.52e-6) )

              DU_bin4 = (LogNormFracBelow(dust_c_mmd,DDspec(DDdep_dustc)%sigma,20.00e-6)&
                       - LogNormFracBelow(dust_c_mmd,DDspec(DDdep_dustc)%sigma,5.00e-6) )

              if(masterproc) &
                write(*,*) 'CloudJ seasalt bin mass fractions: ',SS_bin1, SS_bin2, SS_bin3, SS_bin4, &
                           'CloudJ dust bin mass fractions: ',   DU_bin1, DU_bin2, DU_bin3, DU_bin4
            endif

            AER(L,6)=DUST_F(L) * DU_bin1 ! aerosol [g/m2]
            NAA(L,6)=-6                  ! aerosol index

            AER(L,7)=DUST_C(L) * DU_bin2 ! aerosol [g/m2]
            NAA(L,7)=-7                  ! aerosol index

            AER(L,8)=DUST_C(L) * DU_bin3 ! aerosol [g/m2]
            NAA(L,8)=-8                  ! aerosol index

            AER(L,9)=DUST_C(L) * DU_bin4 ! aerosol [g/m2]
            NAA(L,9)=-9                  ! aerosol index

            AER(L,10)=Forest_F(L) ! aerosol [g/m2]

            ! convert ratio to index matching aerosol carbon-ratio input file
            ! index -10 to -17 fossil fuel, index -18 to -33 for biomas burn
            FF_BC_index = -18 ! first index for biomass burning in UMa.dat

            if(ffirebc_i>0 .and. Forest_F(L)>1e-18) then ! avoid dividing by zero 
              FF_BCratio_real = xn_adv(ffirebc_i,i_emep,j_emep,k) &
                              * 1000 * roa(i_emep,j_emep,k,1) / ATWAIR  &
                              * DeltZ(L) * species_adv(ffirebc_i)%molwt / Forest_F(L)
          
              ! tabulated in 2% increments
              FF_BC_index = max( 0,nint( FF_BCratio_real / 0.02) ) ! round to nearest int
              FF_BC_index = min( FF_BC_index, 15 )
              FF_BC_index = -FF_BC_index - 18 ! from -18 to -33
            endif
            NAA(L,10)=FF_BC_index

          endif ! if uses%cloudjaerosol
    end do ! EMEP temperature, relative humidity, cloud water path & aerosol done

    L=0
    do k=KMAX_BND,1,-1 ! EMEP level interface pressure
          L=L+1
          PPP(L) = A_bnd(k)/100.0 + B_bnd(k)*PSURF ! PSURF in hPa
    end do
    
    L=0 ! EMEP ozone path (molec/cm^2), derived using interface pressure
    do k=KMAX_MID,1,-1
          L=L+1
          OOO(L) = xn_adv(IXADV_O3,i_emep,j_emep,k)*(PPP(L)-PPP(L+1))*MASFAC
    end do

!------------------------------------------------------------------------------------
! Set layers above the EMEP model
!------------------------------------------------------------------------------------

    K = 0 ! populate stratospheric temperature and pressure
    do L=KMAX_BND,OZ_TOP ! KMAX_BND = KMAX_MID + 1
          PPP(L+1) = hpa_ozone(K + satellite_altind)            ! interface
          TTT(L) = temp_obs(i_emep,j_emep,K + satellite_altind) ! mid-lvl
          K = K + 1
    end do

    ! calculate density and geometric altitude (cm) for all levels
    ZZZ(1) = 16.d5*log10(1013.25d0/PPP(1))

    do L = 1,OZ_TOP 
          DDD(L) = (PPP(L)-PPP(L+1))*MASFAC
          SCALEH = 1.3806d-19*MASFAC*TTT(L)
          ZZZ(L+1) = ZZZ(L) -( LOG(PPP(L+1)/PPP(L)) * SCALEH ) 
    end do
   
    K = 0 ! populate stratospheric ozone; convert mol/m3 to molec/cm2
    do L=KMAX_BND,OZ_TOP
          OOO(L) = ozon_obs(i_emep,j_emep,K + satellite_altind) * molec_cm3 * &
                   ( ZZZ(L+1) - ZZZ(L) ) ! lvl thickness
          K = K + 1
    end do
    
    !-----------------------------------------------------------------------
    !---fast-JX:  SOLAR_JX is called only once per grid-square to set U0, etc. 
    call SOLAR_JX(GMTAU,IDAY,YGRD,XGRD,SZA,U0,SOLF)

    LWP(:)      = 0.d0   ! liquid water path (g/m2)
    IWP(:)      = 0.d0   ! ice water path (g/m2)
    AERSP(:,:)  = 0.d0   ! aerosol path (g/m2)
    NDXAER(:,:) = 0      ! aerosol index type
    CLDIW(:)    = 0

    REFFL(:) = 0.d0      ! R-effective(microns) in liquid cloud
    REFFI(:) = 0.d0      ! R-effective(microns) in ice cloud
    
    Dobson(i_emep,j_emep) = 0. ! array to store total dobson units of atmospheric column

    do L = 1,OZ_TOP           
          ! CLDIW is an integer flag: 1 = water cld, 2 = ice cloud, 3 = both
          if (CWLC(L) > 1.d-11) CLDIW(L) = 1
          if (CWIC(L) > 1.d-11) CLDIW(L) = CLDIW(L) + 2

          if (CWLC(L) > 1.d-12) then            ! [kg/kg]
                LWP(L) = CLWP(L)                   ! [g/m2]
                PMID = 0.5d0*(PPP(L)+PPP(L+1))
                F1   = 0.005d0 * (PMID - 610.d0)
                F1   = min(1.d0, max(0.d0, F1))
                REFFL(L) = 9.6d0*F1 + 12.68d0*(1.d0-F1)    
          endif

          if (CWIC(L) > 1.d-12) then            ! [kg/kg]
                IWP(L) = CIWP(L)                   ! [g/m2]
                ZDEL = (ZZZ(L+1) - ZZZ(L))*0.01d0  ! [m]
                ICWC = IWP(L)/ZDEL                 ! [g/m3]
                REFFI(L) = 164.d0 * (ICWC**0.23d0)
          endif

          NDXAER(:,:) = NAA(:,:) ! populate aerosol arrays
          AERSP(:,:)  = AER(:,:)
          
          dobson(i_emep,j_emep) = dobson(i_emep,j_emep) + OOO(L)/2.687e16 ! 1 DU = 2.687e16 molecules of O3 per square centimetre
    end do

    ! albedo read in as %. Get land-name, then read albedo and convert to fraction 
    ALBEDO = 0.0
    do ilu=1, LandCover(i_emep,j_emep)%ncodes
          lu = LandCover(i_emep,j_emep)%codes(ilu)
          ALBEDO = ALBEDO + LandDefs(lu)%Albedo*0.01*LandCover(i_emep,j_emep)%fraction(ilu)
    end do
    REFLB = ALBEDO
                      
    !--CLOUD_JX:   different cloud schemes (4:8 require max-ran overlap algorithm)
    !       CLDFLAG = 1  :  Clear sky J's
    !       CLDFLAG = 2  :  Averaged cloud cover
    !       CLDFLAG = 3  :  cloud-fract**3/2, then average cloud cover
    !       CLDFLAG = 4  :  Average direct solar beam over all ICAs, invert to get clouds
    !       CLDFLAG = 5  :  Random select NRANDO ICA's from all(Independent Column Atmos.)
    !       CLDFLAG = 6  :  Use all (up to 4) quadrature cloud cover QCAs (mid-pts of bin)
    !       CLDFLAG = 7  :  Use all (up to 4) QCAs (average clouds within each Q-bin) ! recommended Prather 2015
    !       CLDFLAG = 8  :  Calcluate J's for ALL ICAs (up to 20,000 per cell!)
    CLDFLAG = 3        
    
    ! asymmetry factor equivalent isotropy: CLDFLAG = 3 only   
    FG0     = 1.1

    ! for CLDFLAG > 3 only
    CLDCOR  = 0.33
    NRANDO  = 5
    LNRG    = 6

    if (first_call) LPRTJ = USES%CLOUDJVERBOSE
    
    !=======================================================================
    ! outputs Jvalues (VALJXX), NICA and JCOUNT
    !
    ! inputs PPP in hPa, ZZZ in cm, TTT in Kelvin, DDD in molec/cm2,
    ! RRR as RH (0.00->1.00), OOO in molec/cm2, LWP in g/m2, IWP in g/m2
    !=======================================================================
    call CLOUD_JX (U0,SZA,REFLB,SOLF,FG0,LPRTJ,PPP,ZZZ,TTT,            &
                   DDD,RRR,OOO,LWP,IWP,REFFL,REFFI,CLF,CLDCOR,CLDIW,   &
                   AERSP,NDXAER,OZ_TOP,naerosol,VALJXX,JVN_,           &
                   CLDFLAG,NRANDO,IRAN,LNRG,NICA,JCOUNT)
    !=======================================================================
        
    ! populate 3D Jvalue array with CloudJ values for the photochemical reactions in EMEP    
    do L=1,KMAX_BND-KCHEMTOP 
          ! NRATJ is number of phot reactions from j2j input file
          do J=1,NRATJ
                if(JIND(J)>0) rcphotslice(J,kmax_bnd-L,i_emep,j_emep) = VALJXX(L,JIND(J))*JFACTA(J)
          end do
    end do

    if (first_call) then 
      do I=1,NRATJ
        ! name to search for must be exactly 10 characters long; first ten characters in FJX_j2j.dat
        if ('O3        ' .eq. JLABEL(I)(:10) ) IDO3_O3P  = I
        if ('O3(1D)    ' .eq. JLABEL(I)(:10) ) IDO3_O1D  = I
        if ('NO2       ' .eq. JLABEL(I)(:10) ) IDNO2     = I   
        if ('H2COa     ' .eq. JLABEL(I)(:10) ) IDHCHO_H  = I
        if ('H2COb     ' .eq. JLABEL(I)(:10) ) IDHCHO_H2 = I 
        if ('H2O2      ' .eq. JLABEL(I)(:10) ) IDH2O2    = I   
        if ('CH3OOH    ' .eq. JLABEL(I)(:10) ) IDCH3O2H  = I   
        if ('NO3c      ' .eq. JLABEL(I)(:10) ) IDNO3     = I ! lumped EMEP
        if ('HNO4      ' .eq. JLABEL(I)(:10) ) IDHO2NO2  = I
        if ('CH3COCHO  ' .eq. JLABEL(I)(:10) ) IDRCOCHO  = I
        if ('BIACET    ' .eq. JLABEL(I)(:10) ) IDCH3COY  = I
        if ('CH3COC2H5 ' .eq. JLABEL(I)(:10) ) IDMEK     = I
        if ('CH3CHO    ' .eq. JLABEL(I)(:10) ) IDCH3CHO  = I
        if ('CHOCHOa   ' .eq. JLABEL(I)(:10) ) IDGLYOXA  = I
        if ('CHOCHOb   ' .eq. JLABEL(I)(:10) ) IDGLYOXB  = I
        if ('CHOCHOc   ' .eq. JLABEL(I)(:10) ) IDGLYOXC  = I
        if ('PAN       ' .eq. JLABEL(I)(:10) ) IDPAN     = I ! only in CJX
        if ('HNO3      ' .eq. JLABEL(I)(:10) ) IDHNO3    = I
        if ('HNO2      ' .eq. JLABEL(I)(:10) ) IDHONO    = I
        if ('GLYOX     ' .eq. JLABEL(I)(:10) ) IDCHOCHO  = I ! lumped EMEP 
        if ('NO3a      ' .eq. JLABEL(I)(:10) ) IDNO3_NO  = I 
        if ('NO3b      ' .eq. JLABEL(I)(:10) ) IDNO3_NO2 = I 
        if ('CH3COCH3a ' .eq. JLABEL(I)(:10) ) IDACETON  = I ! only a-channel
        if ('N2O5      ' .eq. JLABEL(I)(:10) ) IDN2O5    = I

        ! duplicates with different names for historical reasons
        if ('CHOCHOa   ' .eq. JLABEL(I)(:10) ) IDCHOCHO_2CHO = I
        if ('CHOCHOb   ' .eq. JLABEL(I)(:10) ) IDCHOCHO_2CO  = I
        if ('CHOCHOc   ' .eq. JLABEL(I)(:10) ) IDCHOCHO_HCHO = I

        ! non-EmChem photolysis rates
        if ('CH3COCH3a ' .eq. JLABEL(I)(:10) ) IDCH3COCH3  = I
        if ('MCM15     ' .eq. JLABEL(I)(:10) ) MCM_J15     = I
        if ('MCM17     ' .eq. JLABEL(I)(:10) ) MCM_J17     = I
        if ('MeAcr     ' .eq. JLABEL(I)(:10) ) MCM_J18     = I 
        if ('HPALD1    ' .eq. JLABEL(I)(:10) ) MCM_J20     = I
        if ('CH3COC2H5 ' .eq. JLABEL(I)(:10) ) MCM_J22     = I
        if ('MVK       ' .eq. JLABEL(I)(:10) ) MCM_J23     = I
        if ('IPRNO3    ' .eq. JLABEL(I)(:10) ) IDiC3H7ONO2 = I 
        if ('C2H5CHO   ' .eq. JLABEL(I)(:10) ) IDC2H5CHO   = I
        if ('HAC       ' .eq. JLABEL(I)(:10) ) IDACETOL    = I
        if ('GLYC      ' .eq. JLABEL(I)(:10) ) IDGLYALD    = I
        if ('MCRENOL   ' .eq. JLABEL(I)(:10) ) IDMCRENOL   = I
        if ('ETP       ' .eq. JLABEL(I)(:10) ) IDETP       = I
        if ('ETHP      ' .eq. JLABEL(I)(:10) ) IDETHP      = I
        if ('ATOOH     ' .eq. JLABEL(I)(:10) ) IDATOOH     = I
        if ('R4P       ' .eq. JLABEL(I)(:10) ) IDR4P       = I
        if ('RIPC      ' .eq. JLABEL(I)(:10) ) IDRIPC      = I
        if ('PRALDP    ' .eq. JLABEL(I)(:10) ) IDPRALDP    = I
        if ('IDHPE     ' .eq. JLABEL(I)(:10) ) IDIDHPE     = I
        if ('PIP       ' .eq. JLABEL(I)(:10) ) IDPIP       = I
        if ('ITCNa     ' .eq. JLABEL(I)(:10) ) IDITCN      = I
        if ('INPDa     ' .eq. JLABEL(I)(:10) ) IDINPD      = I
        if ('MAP       ' .eq. JLABEL(I)(:10) ) IDMAP       = I
        if ('RP        ' .eq. JLABEL(I)(:10) ) IDRP        = I
        if ('MENO3     ' .eq. JLABEL(I)(:10) ) MCM_J51     = I
        if ('ETNO3     ' .eq. JLABEL(I)(:10) ) MCM_J52     = I
        if ('NPRNO3    ' .eq. JLABEL(I)(:10) ) MCM_J53     = I
        if ('IPRNO3    ' .eq. JLABEL(I)(:10) ) MCM_J54     = I
        if ('PROPNN    ' .eq. JLABEL(I)(:10) ) MCM_J56     = I
        if ('R4N2      ' .eq. JLABEL(I)(:10) ) IDR4N2      = I
        if ('MVKN      ' .eq. JLABEL(I)(:10) ) IDMVKN      = I
        if ('INPB      ' .eq. JLABEL(I)(:10) ) IDINPB      = I
        if ('IHN3      ' .eq. JLABEL(I)(:10) ) IDIHN3      = I
      enddo ! NRATJ

      ! place indices in the photol_used array
      call setPhotolUsed()
  
      ! verify that all phot rates were found (i.e. greater than 0)
      do I=1,NPHOTOLRATES
        if ( photol_used(I) < 0  .and. MasterProc ) then 
          write(*,*) 'Missing photolysis index in photol_used:', I
          call StopAll( 'Photolysis reaction not found in CloudJ input file.')
        endif
      enddo
    endif

    first_call=.false.
    LPRTJ=.false.

end subroutine setup_phot_cloudj


subroutine write_jvals(i_emep,j_emep)

    ! this routine writes the requested 3D J-value output. Indices are set by 
    ! setup_phot_cloudj (CloudJ) 
    !
    ! there is surely a prettier way to do this....but this also works...

    integer, intent(in) :: i_emep,j_emep
    
    ! CloudJ (_fj) photolysis rate output

    if(photo_out_ix_no2_fj>0)then
          d_3d(photo_out_ix_no2_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
            rcphot(IDNO2,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_o3a_fj>0)then
          d_3d(photo_out_ix_o3a_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
            rcphot(IDO3_O3P,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_o3b_fj>0)then
          d_3d(photo_out_ix_o3b_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
            rcphot(IDO3_O1D,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_h2o2_fj>0)then
      d_3d(photo_out_ix_h2o2_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
        rcphot(IDH2O2,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_hno3_fj>0)then
      d_3d(photo_out_ix_hno3_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
        rcphot(IDHNO3,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_ach2o_fj>0)then
      d_3d(photo_out_ix_ach2o_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
        rcphot(IDHCHO_H,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_bch2o_fj>0)then
      d_3d(photo_out_ix_bch2o_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
        rcphot(IDHCHO_H2,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_hono_fj>0)then
      d_3d(photo_out_ix_hono_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
        rcphot(IDHONO,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_ho2no2_fj>0)then
      d_3d(photo_out_ix_ho2no2_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
        rcphot(IDHO2NO2,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_no3_fj>0)then
      d_3d(photo_out_ix_no3_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
        rcphot(IDNO3,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_ch3o2h_fj>0)then
      d_3d(photo_out_ix_ch3o2h_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
        rcphot(IDCH3O2H,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_MEK_fj>0)then
      d_3d(photo_out_ix_MEK_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
        rcphot(IDMEK,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_N2O5_fj>0)then
      d_3d(photo_out_ix_N2O5_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
        rcphot(IDN2O5,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_GLYOX_fj>0)then
      d_3d(photo_out_ix_GLYOX_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
        rcphot(IDCHOCHO,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_CH3CHO_fj>0)then
      d_3d(photo_out_ix_CH3CHO_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
        rcphot(IDCH3CHO,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_ACETON_fj>0)then
      d_3d(photo_out_ix_ACETON_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
        rcphot(IDACETON,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_MGLYOX_fj>0)then
      d_3d(photo_out_ix_MGLYOX_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
        rcphot(IDRCOCHO,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_BIACET_fj>0)then
      d_3d(photo_out_ix_BIACET_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
        rcphot(IDCH3COY,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_PAN_fj>0)then
      d_3d(photo_out_ix_PAN_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
        rcphot(IDPAN,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_GLYOXA_fj>0)then
      d_3d(photo_out_ix_GLYOXA_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
        rcphot(IDGLYOXA,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_GLYOXB_fj>0)then
      d_3d(photo_out_ix_GLYOXB_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
        rcphot(IDGLYOXB,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

    if(photo_out_ix_GLYOXC_fj>0)then
      d_3d(photo_out_ix_GLYOXC_fj,i_emep,j_emep,1:num_lev3d,IOU_INST) = &
        rcphot(IDGLYOXC,max(KCHEMTOP,lev3d(1:num_lev3d))) 
    endif

end subroutine write_jvals


endmodule CloudJ_mod














