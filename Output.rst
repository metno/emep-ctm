.. _`ch-output`:

Output files
============

Output files from a model run are written out in either ASCII,
or (for most data outputs) in NetCDF format.
The different NetCDF files are named after the ``runlabel1`` parameter set in ``modrun.sh``.
The model output is written to the same directory as where the runscript
where submitted, as described in :numref:`ch-submitarun`.

To check your model run, already prepared model result files can be
downloaded using the catalog tool (:numref:`sec-ModelCode`) as follows:

.. code-block:: bash

    # download the output
    catalog.py -R 5.8 --output

Unpacked files are placed in an output directory with model run results
for a whole year, and sometimes with a smaller test run for e.g. April.


.. csv-table:: List of model output files
   :name: tab-outputs
   :header: **Output data files**, **Short description**, **Format**
   :delim: &

   ``Base_hourInst.nc`` & Gridded hourly instantaneous values of a & NetCDF
                        & selection of compounds               &
   ``Base_hour.nc``     & Gridded hourly mean values of a      & NetCDF
                        & selection of compounds               &
   ``Base_day.nc``      & Gridded daily   values of a          & NetCDF
                        & selection of compounds               &
   ``Base_month.nc``    & Gridded monthly values of a          & NetCDF
                        & selection of compounds               &
   ``Base_fullrun.nc``  & Gridded yearly  values of a          & NetCDF
                        & selection of compounds               &
   ``sites_YYYY.nc``    & Single-level station samples, means or totals & NetCDF [#YY]_
   ``sondes_YYYY.nc``   & Instantaneous vertical profiles              & NetCDF [#YY]_
   ``sites_YYYY.csv``   & Optional station text output                 & CSV
   ``sondes_YYYY.csv``  & Optional flattened profile text output       & CSV
   **Additional files**
   ``RunLog.out``       & Summary log of runs, including total emissions & ASCII
                        & of different air pollutants per country        &
   ``Timing_Baser.out`` & Timing log file                                & ASCII
   ``MassBudgetSummary.txt`` & Mass balance summary                      & ASCII

.. rubric:: Footnotes
.. [#YY] ``YYYY``: year.
.. [#Old] Deprecated output.


.. _`sec-outputparam`:

Output parameters NetCDF files
------------------------------

Parameters to be written out into ``Base_hour.nc``, ``Base_day.nc``, ``Base_month.nc`` and
``Base_year.nc`` are defined in ``My_Derived_mod.f90`` and ``Derived_ml.f90``.
In ``My_Derived_mod.f90``, the use can specify the output species (air
concentrations, depositions, column values), units and temporal
resolution of the outputs (daily, monthly, yearly).

The name of output parameter provides some information about data. The
names start with TYPE of the parameter, namely SURF (surface air
concentrations), DDEP (Dry deposition), WDEP (Wet deposition), COLUMN
(Vertically integrated parameters), Area (Surface area) etc.

For surface air concentrations, the general name pattern is
``SURF_UNITS_COMPONENT``. Here, ``UNITS`` can e.g. be "ug" (|ugm3|\ ),
"ugS" (|ugSm3|\ ), "ugN" (|ugNm3|\ ), or "ppb".
Note that the components are classified either as "SPEC" (species) or "GROUP".
The content of complex GROUP components can be found in ``CM_ChemGroups_mod.f90``.

For column integrated parameters, the names are `COLUMN_COMPONENT_kNLAYERS`,
where `NLAYERS` is the number of layers from model top included in the integration.
The units for column outputs are "ugm2" (\ :math:`\mu g/m^2`\ ),
"mcm2" (\ :math:`molec/m^2`\ ) or "e15mcm2" (\ :math:`10^{15} molec/m^2`\ ).

For dry depositions, given per :math:`1 m^2` of specified landuse,
the names look like ``DDEP_COMPONENT_m2LANDUSE``,
where ``LANDUSE`` can be either a specific landuse type or a cell average.
For wet depositions, the names are ``WDEP_COMPONENT``.
The units for dry and wet depositions are |mgm2|\ , |mgSm2| or |mgNm2|\.

``VG_COMPONENT_LANDUSE`` are the dry deposition velocities on various
landuse types, typically in :math:`cm/s`.

:numref:`tab-outpar` lists most of output parameters, providing additional
explanation to the complex components. For a complete suit of currently
selected output parameters, see provided output NetCDF files, or
``My_Derived_mod.f90`` module.

.. csv-table:: List of output parameters
    :name: tab-outpar
    :header: **Parameter name**, **Short description**, **Comments**
    :delim: &

    **Surface Concentrations**
    ``SURF_ppb_O3``         & |O3|   [ppb]                        & Also available in |ugm3|
    ``SURF_ugN_NO``         & NO     [|ugNm3|\ ]                  & Also available in |ugm3| and ppb
    ``SURF_ugN_NO2``        & |NO2|  [|ugNm3|\ ]                  & Also available in |ugm3| and ppb
    ``SURF_ugN_HNO3``       & |HNO3| [|ugNm3|\ ]                  & Also available in |ugm3| and ppb
    ``SURF_ugN_NH3``        & |NH3|  [|ugNm3|\ ]                  & Also available in |ugm3| and ppb
    ``SURF_ugS_SO2``        & |SO2|  [|ugSm3|\ ]                  & Also available in |ugm3| and ppb
    ``SURF_ug_SO4``         & |SO4|  [|ugm3|\ ]                   &
    ``SURF_ug_NO3_F``       & |NO3|  fine aerosol   [|ugm3|\ ]    & As ammonium nitrate
    ``SURF_ug_NO3_C``       & |NO3|  coarse aerosol [|ugm3|\ ]    & Associated with sea salt and mineral dust
    ``SURF_ug_PNO3``        & |NO3|  total          [|ugm3|\ ]    & Sum of fine and coarse nitrate
    ``SURF_ug_NH4_F``       & |NH4|  fine aerosol   [|ugm3|\ ]    & As ammonium sulphate and ammonium nitrate
    ``SURF_ug_SIA``         & SIA [|ugm3|\ ]                      & Secondary Inorganic Aerosol
    ``SURF_ug_ECFINE``      & EC fine   [|ugm3|\ ]                & Elemental carbon
    ``SURF_ug_ECCOARSE``    & EC coarse [|ugm3|\ ]                & Elemental carbon
    ``SURF_ug_PM_OM25``     & OM fine   [|ugm3|\ ]                & Organic Matter fine aerosol, also available in |ugCm3|
    ``SURF_ug_PM_OMCOARSE`` & OM coarse [|ugm3|\ ]                & Organic Matter coarse aerosol, also available in |ugCm3|
    ``SURF_ug_SEASALT_F``   & Sea salt fine aerosol    [|ugm3|\ ] &
    ``SURF_ug_SEASALT_C``   & Sea salt coarse aerosol  [|ugm3|\ ] &
    ``SURF_ug_SS``          & Sea salt                 [|ugm3|\ ] & Sum of fine and coarse sea salt
    ``SURF_ug_DUST_ROAD_F`` & Road dust fine aerosol   [|ugm3|\ ] &
    ``SURF_ug_DUST_ROAD_C`` & Road dust coarse aerosol [|ugm3|\ ] &
    ``SURF_ug_DUST_WB_F``   & Windblown dust fine      [|ugm3|\ ] &
    ``SURF_ug_DUST_WB_C``   & Windblown dust coarse    [|ugm3|\ ] &
    ``SURF_ug_DUST_SAH_F``  & Saharan dust fine        [|ugm3|\ ] & From Boundary conditions
    ``SURF_ug_DUST_SAH_C``  & Saharan dust coarse      [|ugm3|\ ] & From Boundary conditions
    ``SURF_ug_DUST_NAT_F``  & Natural dust fine        [|ugm3|\ ] & Windblown and Saharan
    ``SURF_ug_DUST_NAT_C``  & Natural dust coarse      [|ugm3|\ ] & Windblown and Saharan
    ``SURF_ug_DUST``        & Mineral dust             [|ugm3|\ ] & From all sources
    ``SURF_ug_PM10``        & |PM10| dry [|ugm3|\ ]               &
    ``SURF_ug_PM10_rh50``   & |PM10| wet [|ugm3|\ ]               & |PM10| particle water at 50 %rh
    ``SURF_ug_PM25``        & |PM25| dry [|ugm3|\ ]               & Includes fine PM
    ``SURF_ug_PM25_rh50``   & |PM25| wet [|ugm3|\ ]               & |PM25| particle water at 50 %rh
    ``SURF_ug_PM25X``       & |PM25| dry [|ugm3|\ ]               & Includes fine PM, EC and OM
    ``SURF_ug_PM25X_rh50``  & |PM25|     [|ugm3|\ ]               & As ``PM25X`` + particle water at 50 %rh
    ``SURF_ug_PMFINE``      & Fine PM [|ugm3|\ ]                  & Sum of all fine aerosols
    ``SURF_ug_PPM25``       & Primary |PM25| [|ugm3|\ ]           & Anthropogenic emissions
    ``SURF_ug_PPM_C``       & Primary coarse PM [|ugm3|\ ]        & Anthropogenic emissions
    ``SURF_ug_PM25_FIRE``   & |PM25| from forest fires [|ugm3|\ ] & Sum of BC, OC and rest |PM25|
     **Dry Depositions**
    ``DDEP_SOX_m2Grid``     & Oxidized sulphur  [|mgSm2|\ ]       & For a grid cell landuse area weighted
    ``DDEP_SOX_m2Conif``    & Oxidized sulphur  [|mgSm2|\ ]       & To coniferous forest
    ``DDEP_NOX_m2Grid``     & Oxidized nitrogen [|mgNm2|\ ]       & For a grid cell landuse area weighted
    ``DDEP_NOX_m2Decid``    & Oxidized nitrogen [|mgNm2|\ ]       & To decideous forest
    ``DDEP_RDN_m2Grid``     & Reduced nitrogen  [|mgNm2|\ ]       & For a grid cell landuse area weighted
    ``DDEP_RDN_m2Seminat``  & Reduced nitrogen  [|mgNm2|\ ]       & To semi-natural
    **Wet Depositions**
    ``WDEP_PREC``           & Precipitation     [mm]              &
    ``WDEP_SOX``            & Oxidized sulphur  [|mgSm2|\ ]       &
    ``WDEP_SS``             & Sea salt          [|mgm2|\ ]        &
    **Others** &&
    ``AOD``                 & Aerosol Optical Depth at 550nm      & Experimental
    ``Area_Crops_Frac``     & Area fraction of crops              & Available for several landuses
    ``VG_NO3_F_Grid``       & Dry deposition velocity of fine |NO3| & Grid cell average
    **Meteorological parameters**
    ``USTAR_GRID``          & :math:`U^*` grid averaged           & Available for several landuses
    ``T2m``                 & Temperature at 2m [|degC|\ ]        &
    ``rh2m``                & Fractional relative humidity at 2m  &



Emission outputs
----------------

``Emis_mgm2_XX`` fields in the output, give all emissions used by the model (accumulated over the relevant period). ``Sec_Emis_mgm2_XX`` are "sector emissions" as defined in the sector emission files, i.e. includes only contributions from the files defined in Emis_sourceFiles. ``Sec_Emis_mgm2_XX`` do not include emissions such as volcanoes, forest fires, DMS, lightning, aircraft etc.
Note that for VOC, some emissions from the sector file are set as "unreactive"; those will show under ``Sec_Emis_mgm2_voc`` but not under ``Emis_mgm2_voc``.
For hourly outputs of emissions set

.. code-block:: fortran

  HourlyEmisOut = T,

For daily outputs of emissions set

.. code-block:: fortran

  DailyEmisOut = T,

Detailed emissions by sectors can be obtained with the keyword ``SecEmisOutWanted`` for the wanted sectors. For example adding the lines:

.. code-block:: fortran

  SecEmisOutWanted(2) = T,
  SecEmisOutWanted(7) = T,

will give you the emissions for sector 2 and 7 for all components.

Totals per country and sectors (all), can be obtained in the log with:

.. code-block:: fortran

  SecEmisTotalsWanted = T,

To get emissions partitioned into splitted compounds (up to 18), the value ``EmisSplit_OUT=T`` must be set in ``config_emep.nml``

Single species sector specific emissions can be obtained with for example:

.. code-block:: fortran

    'BENZENE_GNFR_C','SpecSecEmis','BENZENE','-','mg/m2',3, -99, T , 1.0e6,     F, 'YM',


Add your own fields
-------------------

Most standard output can be outputted by adding lines and modifying the parameters in the ``config_emep.nml`` file.

The meteorological fields defined in the ``met`` array in the ``MetFields_mod.f90`` file, can be retrieved by using the 'MET2D' or 'MET3D' keywords. If a 3D array is requested with the 'MET2D' keyword, only the lowest level is written out.

If you want an array that does not fit in any category, or even make your own special field, you can get it in the output using the procedure shown below; this will however require that you write in the code and recompile.
For instance in config_emep.nml OutputMisc define:

``  'J(NO2)'  ,'USET','D3_J(NO2)'  ,'photorate','1/s' ,-99,-99,F,1.0,T,'H',``

- The first column (name) is the name as shown in the output
- The second column (class) must be 'USET'
- The strings of the first and third columns can be chosen freely, but if one of them starts with the two characters 'D3', it will be interpreted as a 3 dimensional field
- The fourth column can be any string
- The fifth column is the unit, as show in the output
- The sixth column (index) is an integer that can be used to characterize internal indices
- The seventh columns should be a negative integer
- The eigth column can be F or T, indicating wether the field must be divided by the time step (dt_advec)
- The ninth column (scale) is a scaling factor
- The tenth column, F or T, indicates if the field must be averaged (T) or accumulated (F)
- The eleventh (last) column indicates the periodicity of the output. 'H'-> every hour, 'YMH'--> every hour, month and at the end of the run (and other combinations are allowed).

In the code you must define the indice of your new ouput. The requested outputs strings are stored in f_2d and f_3d; for instance

.. code-block:: fortran

    photo_out_ix = find_index("D3_J(NO2)", f_3d(:)%subclass)

and the values of the field must be put into the d_2d or d_3d array, using this index, for instance:

.. code-block:: fortran

    if(photo_out_ix>0) d_3d(photo_out_ix,i,j,1:num_lev3d,IOU_INST) = rcphot(IDNO2,lev3d(1:num_lev3d))

(for 2D output, write in d_2d and ommit the vertical index)

.. _sec-sitessondes-output:

Site time series and sonde profiles
----------------------------------

Station output writes selected model data directly at locations specified
in the files described in :ref:`sec-sitessondes-input`.

Output files
~~~~~~~~~~~

The model writes:

* ``sites_YYYY.nc``: single-level station data.
* ``sondes_YYYY.nc``: vertical profiles at selected locations.
* ``sites_YYYY.csv`` and ``sondes_YYYY.csv``: optional text output.

These filenames do not include ``runlabel1``. Use separate run directories
to avoid overwriting station output from another simulation.

The standard output interval is one hour. Site files can contain
instantaneous samples, interval means and interval totals in the same
file. Sonde profiles are instantaneous samples.

Station output is suppressed in source-receptor mode.

NetCDF structure
~~~~~~~~~~~~~~~

In the usual NetCDF/CDL display order, site variables have dimensions
``(time, station)`` and sonde profile variables have dimensions
``(time, lev, station)``.

Sonde levels are ordered from the lowest model layer upwards. All
``KMAX_MID`` layers are included; the gridded-output ``lev3d`` selection
does not restrict sonde profiles.

Station identifiers and model-grid metadata accompany the data.
Geographical location input supplies latitude and longitude; these are
not inferred from ``IJKdown`` input.

Sonde files contain surface pressure ``PS`` in hPa and hybrid-coordinate
coefficients. Midpoint pressure can be reconstructed as:

.. code-block:: text

    pressure_hPa = hyam + hybm * PS

The optional diagnostic ``p_mid`` writes midpoint pressure directly in Pa.

Sampling and temporal statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Raw advected species
    Instantaneous values in ppb. Site output at the lowest model layer
    includes the existing surface correction. Higher-layer site samples
    and sonde profiles use model-layer concentrations without that
    lowest-layer site correction.

Raw short-lived species
    Instantaneous values in molecules/cm3 at the selected site layer
    or throughout the sonde column.

Direct site diagnostics
    ``T2`` is instantaneous two-metre temperature in degrees Celsius.
    ``th`` is instantaneous potential temperature in K at the selected
    site layer.

Two-dimensional derived site fields
    Fields marked for averaging produce interval means. Accumulated
    fields produce interval totals. Timestep scaling and unit conversion
    follow the derived-field definition.

A mapped two-dimensional field is sampled at the station's horizontal
grid cell. Its own definition determines its vertical interpretation;
the station's selected layer does not make it an arbitrary-height field.

For example, ``O3`` selects instantaneous ozone at the chosen site layer,
whereas ``SURF_ppb_O3`` selects an interval-mean derived surface field.

Statistical products requiring dedicated handling, including supported
rejection cases for daily maxima, rolling ozone metrics and fixed-hour
diagnostics, cannot be obtained through the generic site accumulator.

Sonde diagnostics
~~~~~~~~~~~~~~~~

Besides available chemical species, the direct profile diagnostics 
available are:

.. list-table::
   :header-rows: 1

   * - Name
     - Quantity
     - Units
   * - ``PM25``
     - Fine particulate mass plus the implemented coarse contribution
     - ug/m3
   * - ``PMco``
     - Coarse particulate mass
     - ug/m3
   * - ``NOy``
     - Sum over the implemented oxidized-nitrogen group
     - ppb
   * - ``SH``
     - Specific humidity
     - kg/kg
   * - ``RH``
     - Relative humidity, capped at 1
     - fraction
   * - ``roa``
     - Air density
     - kg/m3
   * - ``dz``
     - Layer thickness
     - m
   * - ``z_mid``
     - Model layer-midpoint height
     - m
   * - ``p_mid``
     - Layer-midpoint pressure
     - Pa
   * - ``Kz_m2s``
     - Vertical diffusivity diagnostic
     - m2/s
   * - ``th``
     - Potential temperature
     - K
   * - ``T``
     - Air temperature
     - K
   * - ``U``, ``V``
     - Physical wind components along the model grid directions
     - m/s

``T`` and ``RH`` can be requested independently and in either order.

``U`` and ``V`` include the model mapping-factor correction. They are
not necessarily geographic eastward and northward components; rotated
or projected grids require an additional rotation for that interpretation.

General three-dimensional derived fields are not supported by the
sonde writer.

Timestamps
~~~~~~~~~

CSV and NetCDF use the same midpoint-oriented timestamp:

.. code-block:: text

    timestamp = output boundary - half the station output interval

Hourly output written at 01:00 UTC is therefore labelled 00:30 UTC.
The shift includes date changes across midnight.

These labels align with gridded hourly-mean labels when both outputs
use the same interval.

.. important::

    A midpoint label does not imply averaging.

    Raw species, direct site diagnostics and sonde profiles remain
    instantaneous samples taken at the output boundary. Their sampling
    time is half an interval later than the stored timestamp.

CSV dates use ``DD/MM/YYYY`` and times use ``hh:mm``.

For a continuing run, the final New Year midnight record remains in
the previous year's file. A first-ever station record at New Year
creates its initial output file.

Spin-up and output start
~~~~~~~~~~~~~~~~~~~~~~

Both destinations wait until ``spinup_enddate`` and ``out_startdate``
have been reached.

Derived site samples are not accumulated during spin-up. Site accumulators
are reset at scheduled output boundaries where writing is suppressed.

The first written mean or total can cover a partial interval, particularly
when spin-up ends between output boundaries. Its label still uses half
the nominal output interval; it does not identify the exact midpoint
of a shortened accumulation period.

If derived site output is requested but no valid samples have been
accumulated, the writer stops rather than silently reporting zero.

CSV output
~~~~~~~~~

CSV output is enabled by default. To retain NetCDF only, set:

.. code-block:: fortran

    WRITE_STATION_CSV = .false.,

CSV files contain a descriptive preamble and station metadata before the
time-series header. Readers must locate that header rather than treating
the first line as an ordinary CSV table header.

Site records contain:

.. code-block:: text

    site,code,date,hh,variable_1,variable_2,...

Sonde records flatten profiles by variable, with levels ordered from
the lowest layer upwards within each variable.

.. note::

    The current sonde CSV header lists variable names without expanding
    them into separate level columns. NetCDF is recommended for profile
    analysis.

CSV numeric formatting preserves fewer digits than NetCDF. It should
not be used to test exact numerical equality at full model precision.

Comparison with gridded output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a like-for-like comparison:

#. Select the same grid cell, rather than a bilinearly interpolated value.
#. Match the field definition and units.
#. Match vertical sampling and surface-correction treatment.
#. Match averaging or accumulation intervals.
#. Decode timestamps and compare the intervals they represent.

For gridded hourly-mean surface ozone, select ``SURF_ppb_O3`` at sites,
not the raw species ``O3``.

For instantaneous comparisons, account for the difference between the
stored midpoint label and actual output-boundary sampling time.

Do not apply the historical extra one-hour correction to revised station
timestamps. Rounding timestamps is not a substitute for matching their
meaning.

Exclude partial initial intervals when testing equality of complete
hourly means or totals.

.. _`ch-output-config`:

Output configuration
====================

Surface concentrations, column integrated, wet and dry deposition outputs
are defined by the user in the ``Model_config``namelist, in the configuration file ``config_emep.nml``.
Surface concentrations and column integrated outputs specified in ``OutputConcs`` variable
as shown in :numref:`output-concs` and :numref:`output-concs-alt`.
Dry and wet deposition outputs specified in ``DDEP_ECOS``, ``DDEP_WANTED`` and ``WDEP_WANTED`` variables
as shown in :numref:`output-deps`.

.. code-block:: Fortran
  :name: output-concs
  :caption: Output configuration exmple for surface concentrations correspondig to :numref:`tab-outpar`.

    OutputConcs=
      'O3'                  ,'ppb' ,'2d','AIR_CONCS','SPEC' ,'YMDI',
      'NO'                  ,'ugN' ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'NO2'                 ,'ugN' ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'HNO3'                ,'ugN' ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'NH3'                 ,'ugN' ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'SO2'                 ,'ugS' ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'SO4'                 ,'ug'  ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'NO3_F'               ,'ug'  ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'NO3_C'               ,'ug'  ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'PNO3'                ,'ug'  ,'2d','AIR_CONCS','GROUP','YMD',
      'NH4_F'               ,'ug'  ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'SIA'                 ,'ug'  ,'2d','AIR_CONCS','GROUP','YMD',
      'ECFINE'              ,'ug'  ,'2d','AIR_CONCS','GROUP','YMD',
      'ECCOARSE'            ,'ug'  ,'2d','AIR_CONCS','GROUP','YMD',
      'OM25'               ,'ug_PM','2d','AIR_CONCS','GROUP','YMD',
      'OMCOARSE'           ,'ug_PM','2d','AIR_CONCS','GROUP','YMD',
      'SEASALT_F'           ,'ug'  ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'SEASALT_C'           ,'ug'  ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'SS'                  ,'ug'  ,'2d','AIR_CONCS','GROUP','YMD',
      'DUST_ROAD_F'         ,'ug'  ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'DUST_ROAD_C'         ,'ug'  ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'DUST_WB_F'           ,'ug'  ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'DUST_WB_C'           ,'ug'  ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'DUST_SAH_F'          ,'ug'  ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'DUST_SAH_C'          ,'ug'  ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'DUST_NAT_F'          ,'ug'  ,'2d','AIR_CONCS','GROUP','YMD',
      'DUST_NAT_C'          ,'ug'  ,'2d','AIR_CONCS','GROUP','YMD',
      'DUST'                ,'ug'  ,'2d','AIR_CONCS','GROUP','YMD',
      'SURF_PM25water'      ,'ug'  ,'2d','PM25water','MISC' ,'',
      'PMFINE'              ,'ug'  ,'2d','AIR_CONCS','GROUP','YMD',
      'SURF_ug_PM25'        ,'ug'  ,'2d','PM25'     ,'MISC' ,'YMD',
      'SURF_ug_PM25X'       ,'ug'  ,'2d','PM25X'    ,'MISC' ,'YMD',
      'SURF_ug_PM25X_rh50'  ,'ug'  ,'2d','PM25X_rh50','MISC','YMD',
      'SURF_ug_PM25_rh50'   ,'ug'  ,'2d','PM25_rh50','MISC' ,'YMD',
      'PMCO'                ,'ug'  ,'2d','AIR_CONCS','GROUP','',
      'PM10'                ,'ug'  ,'2d','AIR_CONCS','GROUP','YMD',
      'SURF_ug_PM10_rh50'   ,'ug'  ,'2d','PM10_rh50','MISC' ,'YMD',
      'PPM25'               ,'ug'  ,'2d','AIR_CONCS','GROUP','YMD',
      'PPM_C'               ,'ug'  ,'2d','AIR_CONCS','GROUP','YMD',
      'PPM25_FIRE'          ,'ug'  ,'2d','AIR_CONCS','GROUP','YMD',

.. code-block:: Fortran
  :name: output-concs-alt
  :caption: Output configuration exmple for surface concentrations in alternative units.

    OutputConcs=
      'O3'                  ,'ug'  ,'2d','AIR_CONCS','SPEC' ,'YMDI',
      'NO'                  ,'ppb' ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'NO'                  ,'ug'  ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'NO2'                 ,'ppb' ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'NO2'                 ,'ug'  ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'HNO3'                ,'ug'  ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'HNO3'                ,'ppb' ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'NH3'                 ,'ug'  ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'NH3'                 ,'ppb' ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'SO2'                 ,'ug'  ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'SO2'                 ,'ppb' ,'2d','AIR_CONCS','SPEC' ,'YMD',
      'OM25'              ,'ugC_PM','2d','AIR_CONCS','GROUP','YMD',
      'OMCOARSE'          ,'ugC_PM','2d','AIR_CONCS','GROUP','YMD',

.. code-block:: Fortran
  :name: output-deps
  :caption: Output configuration correspondig to :numref:`tab-outpar` dry and wet depositions.

    DDEP_ECOS =
      'Grid'     ,'YMD',
      'Conif'    ,'YM',
      'Seminat'  ,'YM',
      'Water_D'  ,'Y',
      'Decid'    ,'Y',
      'Crops'    ,'Y',

    DDEP_WANTED =
      'SOX'      ,'GROUP','mgS',
      'OXN'      ,'GROUP','mgN',
      'RDN'      ,'GROUP','mgN',

    WDEP_WANTED =
      'PREC'     ,'PREC' ,'mm' ,'YMD',
      'SOX'      ,'GROUP','mgS','YMD',
      'OXN'      ,'GROUP','mgN','YMD',
      'RDN'      ,'GROUP','mgN','YMD',
      'SO2'      ,'SPEC' ,'mgS','YMD',
      'HNO3'     ,'SPEC' ,'mgN','YMD',
