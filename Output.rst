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
    catalog.py -R 5.0 --output

Unpacked files are placed in an output directory with model run results
for a whole year, and sometimes with a smaller test run for e.g. April.


.. csv-table:: List of model output files
   :name: tab-outputs
   :header: **Output data files**, **Short description**, **Format**
   :delim: &

   ``Base_hour.nc``     & Gridded hourly  values of a          & NetCDF
                        & selection of compounds               &
   ``Base_day.nc``      & Gridded daily   values of a          & NetCDF
                        & selection of compounds               &
   ``Base_month.nc``    & Gridded monthly values of a          & NetCDF
                        & selection of compounds               &
   ``Base_fullrun.nc``  & Gridded yearly  values of a          & NetCDF
                        & selection of compounds               &
   ``sites_YYYY.nc``    & Surface daily   values of a          & NetCDF [#YY]_
                        & selection of stations and compounds  &
   ``sondes_YYYY.nc``   & Vertical daily  values of a          & NetCDF [#YY]_
                        & selection of stations and compounds  &
   ``sites_YYYY.cvs``   & ASCII version of ``sites_YYYY.nc``   & ASCII
   ``sondes_YYYY.csv``  & ASCII version of ``sondes_YYYY.nc``  & ASCII [#Old]_
   **Additional files**
   ``RunLog.out``       & Summary log of runs, including total emissions & ASCII
                        & of different air pollutants per country        &
   ``Timing.out``       & Timing log file                                & ASCII

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

Surface concentrations, column integrated, wet and dry deposition outputs
are defined by the user in the ``Model_config``namelist, in the configuration file ``config_emep.nml``.
Surface concentrations and column integrated outputs specified in ``OutputConcs`` variable.
Dry and wet deposition outputs specified in ``DDEP_ECOS``, ``DDEP_WANTED`` and ``WDEP_WANTED`` variables.

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
    ``SURF_ppb_O3``         & |O3|   [ppb]                        &
    ``SURF_ugN_NO``         & NO     [|ugNm3|\ ]                  & Available also in ppb
    ``SURF_ugN_NO2``        & |NO2|  [|ugNm3|\ ]                  & Available also in ppb
    ``SURF_ugN_HNO3``       & |HNO3| [|ugNm3|\ ]                  & Available also in ppb
    ``SURF_ugN_NH3``        & |NH3|  [|ugNm3|\ ]                  & Available also in ppb
    ``SURF_ugS_SO2``        & |SO2|  [|ugSm3|\ ]                  & Available also in ppb
    ``SURF_ug_SO4``         & |SO4|  [|ugm3|\ ]                   &
    ``SURF_ug_NO3_F``       & |NO3|  fine aerosol   [|ugm3|\ ]    & As ammonium nitrate
    ``SURF_ug_NO3_C``       & |NO3|  coarse aerosol [|ugm3|\ ]    & Associated with sea salt and mineral dust
    ``SURF_ug_TNO3``        & |NO3|  total          [|ugm3|\ ]    & Sum of fine and coarse nitrate
    ``SURF_ug_NH4_F``       & |NH4|  fine aerosol   [|ugm3|\ ]    & As ammonium sulphate and ammonium nitrate
    ``SURF_ug_SIA``         & SIA [|ugm3|\ ]                      & Secondary Inorganic Aerosol
    ``SURF_ug_SIA``         & SIA       [|ugm3|\ ]                & Secondary Inorganic Aerosol
    ``SURF_ug_ECFINE``      & EC fine   [|ugm3|\ ]                & Elemental carbon
    ``SURF_ug_ECCOARSE``    & EC coarse [|ugm3|\ ]                & Elemental carbon
    ``SURF_ug_PM_OM25``     & OM fine   [|ugm3|\ ]                & Organic Matter fine aerosol
    ``SURF_ug_PM_OMCOARSE`` & OM coarse [|ugm3|\ ]                &  Organic Matter coarse aerosol
    ``SURF_ug_SEASALT_F``   & Sea salt fine aerosol    [|ugm3|\ ] &
    ``SURF_ug_SEASALT_C``   & Sea salt coarse aerosol  [|ugm3|\ ] &
    ``SURF_ug_SEASALT``     & Sea salt                 [|ugm3|\ ] & Sum of fine and coarse sea salt
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
    ``SURF_ug_PM25``        & |PM25| dry [|ugm3|\ ]               & Includes fine PM and 27% of coarse |NO3|
    ``SURF_ug_PM25_rh50``   & |PM25| wet [|ugm3|\ ]               & |PM25| particle water at 50 %rh
    ``SURF_ug_PM25X``       & |PM25| dry [|ugm3|\ ]               & Includes fine PM and 27% of coarse |NO3|\ , EC and OM
    ``SURF_ug_PM25X_rh50``  & |PM25|     [|ugm3|\ ]               & As ``PM25X`` + particle water at 50 %rh
    ``SURF_ug_PMFINE``      & Fine PM [|ugm3|\ ]                  & Sum of all fine aerosols
    ``SURF_ug_PPM25``       & Primary P|PM25| [|ugm3|\ ]          & Anthropogenic emissions
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

``Emis_mgm2_XX`` fields in the output, give all emissions used by the model (accumulated over the relevant period). ``Sec_Emis_mgm2_XX`` are "sector emissions", i.e. includes only contributions from the files defined in emis_inputlist and Emis_sourceFiles. ``Sec_Emis_mgm2_XX`` do not include emissions such as volcanoes, forest fires, DMS, lightning, aircraft etc.

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

To get emissions partitioned into splitted compounds (up to 18), the value ``EmisSplit_OUT=.true.`` must be set in ``Config_module.f90``, and the code recompiled. (This parameter cannot be set in ``config_emep.nml`` for now)


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

.. _`sec-sitessondes-output`:

ASCII outputs: sites and sondes
-------------------------------

As noted in :numref:`ec-sitessondes-input`, two
main options are available for the output of ASCII files for
comparison with measurements or detailed model analysis at specific sites. These are

sites
    output of surface concentrations for a set of specified measurement
    site locations.

sondes
    output of concentrations for the vertical column above a set of
    specified locations.

The output files ``sites_2015.csv`` and ``sondes_2015.csv`` are comma
separated files that can be read by excel, python or fortran tools.
If you include the whole year, or the 31\ :sup:`st` December,
``sites_2016.csv`` and ``sondes_2016.csv`` are also included in the output.
