!-------------------------------------------------------------------------!
!     Mapping of emissions from GFASv1 to EMEP EmChem16z species
!     Included by  ForestFire_ml
!-------------------------------------------------------------------------!

  character(len=*),  parameter :: &
    BiomassBurningMapping = "GFASv1toEmChem16z"
  integer,  parameter :: &
    NBB_DEFS   = 25-2, & ! No mapping lines below
    NEMEPSPECS = 20-1    ! No EMEP chemical mech specs used

!-----------------------------------------------------------------------!
! Column-1 gives GFASv1 species,
! Column-2 gives MW if needed to get to kg/day. Not needed for GFASv1
!  (GFAS units are [kg/m2/s]. No unit conversion is needed.)
! Column-3 gives mass fraction of GFASv1 emissions assign to EMEP species
! Column-4 gives EMEP species
!-----------------------------------------------------------------------!
  type(bbtype),  dimension(NBB_DEFS) :: FF_defs = (/ &
    bbtype("cofire       ", 1.000, 1.000, CO            ), & ! 081.210 | Carbon Monoxide
    bbtype("ch4fire      ", 1.000, 1.000, CH4           ), & ! 082.210 | Methane
    bbtype("h2fire       ", 1.000, 1.000, H2            ), & ! 084.210 | Hydrogen
    bbtype("noxfire      ", 1.000, 1.000, NO            ), & ! 085.210 | Nitrogen Oxide. Use as NO
!   bbtype("pm10fire     ", 1.000, 1.000, ffire_rempPMco), & ! 087.210 | PM 10
!   bbtype("pm2p5fire    ", 1.000,-1.000, ffire_rempPMco), & ! 087.210 | PM 10 - PM 2.5
    bbtype("pm2p5fire    ", 1.000, 1.000, ffire_rempPM25), & ! 087.210 | PM 2.5
    bbtype("ocfire       ", 1.700,-1.000, ffire_rempPM25), & ! 090.210 | PM 2.5 - Organic Matter
    bbtype("bcfire       ", 1.000,-1.000, ffire_rempPM25), & ! 091.210 | PM 2.5 - Black Carbon
    bbtype("ocfire       ", 1.700, 1.000, ffire_OM      ), & ! 090.210 | Organic Carbon --> O.Matter (OM/OC=1.7)
    bbtype("bcfire       ", 1.000, 1.000, ffire_BC      ), & ! 091.210 | Black Carbon
    bbtype("so2fire      ", 1.000, 1.000, SO2           ), & ! 102.210 | Sulfur Dioxide
    bbtype("ch3ohfire    ", 1.000, 1.000, CH3OH         ), & ! 103.210 | Methanol
    bbtype("c2h5ohfire   ", 1.000, 1.000, C2H5OH        ), & ! 104.210 | Ethanol
    bbtype("c2h4fire     ", 1.000, 1.000, C2H4          ), & ! 106.210 | Ethene
    bbtype("c3h6fire     ", 1.000, 1.000, C3H6          ), & ! 107.210 | Propene
    bbtype("c5h8fire     ", 1.000, 1.000, C5H8          ), & ! 108.210 | Isoprene
    bbtype("toluenefire  ", 1.000, 1.000, TOLUENE       ), & ! 110.210 | Toluene lump(C7H8+C6H6+C8H10)
    bbtype("hialkenesfire", 1.000, 1.000, C3H6          ), & ! 111.210 | Higher Alkenes (CnH2n, C>=4)
    bbtype("hialkanesfire", 1.000, 1.000, NC4H10        ), & ! 112.210 | Higher Alkanes (CnH2n+2, C>=4)
    bbtype("ch2ofire     ", 1.000, 1.000, HCHO          ), & ! 113.210 | Formaldehyde
    bbtype("c2h4ofire    ", 1.000, 1.000, CH3CHO        ), & ! 114.210 | Acetaldehyde
    bbtype("nh3fire      ", 1.000, 1.000, NH3           ), & ! 116.210 | Ammonia
    bbtype("c2h6fire     ", 1.000, 1.000, C2H6          ), & ! 118.210 | Ethane
    bbtype("c4h10fire    ", 1.000, 1.000, NC4H10        ) /) ! 238.210 | Butanes
