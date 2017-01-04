.. _`ch-output`:

Output filesOutput files
========================

Output files from a model run are written out in either ASCII,
or (for most data outputs) in NetCDF format.
The different NetCDF files are named after the ``runlabel1`` parameter set in ``modrun.sh``.
The model output is written to the same directory as where the runscript
where submitted, as described in :numref:`ch-submitarun`.

To check your model run, already prepared model result files can be
downloaded using the `catalog tool`_ (see :numref:`sec-ModelCode`) as follows:

.. code-block:: bash

    # download the output
    catalog.py -R rv4_10 --output

Unpacked files are placed in an output
directory with model run results for a whole year and sometimes with a
smaller test run for i.e. April.

.. _`tab-outputs`:

.. csv-table:: List of model output files
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
   ``sites_YYYY.cvs``   & ASCII version of ``sites_YYYY.nc``   & ASCII [#Old]_
   ``sondes_YYYY.csv``  & ASCII version of ``sondes_YYYY.nc``  & ASCII [#Old]_
   **Additional files**
   ``RunLog.out``       & Summary log of runs, including total emissions & ASCII
                        & of different air pollutants per country        &
    ``Timing.out``      & Timing log file                                & ASCII

.. rubric:: Footnotes
.. [#YY] ``YYYY``: year.
.. [#Old] Deprecated output.


.. _`sec-outputparam`:

Output parameters NetCDF files
------------------------------

Parameters to be written out ``Base_day.nc``, ``Base_month.nc`` and
``Base_year.nc`` are defined in ``My_Derived_ml.f90`` and ``Derived_ml.f90``.
In ``My_Derived_ml.f90``, the use can specify the output species (air
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
The content of complex GROUP components can be found in ``CM_ChemGroups_ml.f90``.

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
are defined by the user in ``config_emep.nml`` file.
Surface concentrations and column integrated outputs
are described in ``OutputConcs_config`` namelist,
Dry and wet deposition outputs
are described in ``OutputDep_config`` namelist.

``VG_COMPONENT_LANDUSE`` are the dry deposition velocities on various
landuse types, typically in :math:`cm/s`.

:numref:`tab-outpar` lists most of output parameters, providing additional
explanation to the complex components. For a complete suit of currently
selected output parameters, see provided output NetCDF files, or
``My_Derived_ml.f90`` module.

.. _`tab-outpar`:

.. csv-table:: List of output parameters
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


.. _`sec-sitesonde`:

ASCII outputs: sites and sondes
-------------------------------

Two main options are available for the output of ASCII files for
comparison with measurements or detailed model analysis. These are

sites
    output of surface concentrations for a set of specified measurement
    site locations.

sondes
    output of concentrations for the vertical column above a set of
    specified locations.

Both sites and sondes are specified and handled in similar ways, in the
module ``Sites_ml.f90``, so we treat them both together below.
Locations are specified in the input files ``sites.dat`` and ``sondes.dat``.
The files start with a description of its content followed by a list of
the stations. For example, a sondes.dat input file may look like this:

.. literalinclude:: sites.dat
    :caption: Minimum ``modrun.sh`` example.

The first line in each file is a header with file content. Then, the
contents are described in more detail. Text strings after ``#`` are just
clarifying comments. 'Area', e.g., is the domain to which the stations
belong, e.g. 'Northern Hemisphere'.

Text after ``:`` is read in by the model:

Units
    Either 'deg' (degrees) or 'index' (model grid indices).

Coords
    Either 'LatLong' (latitudes/longitudes)
    or 'ModelCoords' (indices of the grid box in which the station is located).

VertCoords
    Vertical coordinate system that is used in the model (usually 'EMEPsigma').

Both ``sites.dat`` and ``sondes.dat`` files are optional, but recommended.
The species and meteorological data requested for site and sone output are
specified in ``My_Outputs.f90`` by the use of arrays.
Only a few met fields are defined so far but more can be added into
``Sites_ml.f90`` as required.

The output files ``sites_2014.csv`` and ``sondes_2010.csv`` are comma
separated files that can be read by excel.
If you include the whole year, or the 31\ :sup:`st` December,
``sites_2015.csv`` and ``sondes_2015.csv`` are also incued in the output.
