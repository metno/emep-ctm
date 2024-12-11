.. _`ch-inputfiles`:

Input files
===========

This chapter provides an overview on the necessary input files to run
the EMEP/MSC-W model. A complete set of input files is provided as part of the
EMEP/MSC-W Open Source release to allow model runs for the meteorological
year 2015. :numref:`tab-inputdata` lists the input files.

In the latest release, meteorology is provided for 2 different model domains and resolutions:

- `EECCA` domain with a horizontal resolution of 50x50 km2 (at 60°N), 
  on polar stereographic projection, and 20 vertical levels;
- `EMEP01` domain with a 0.1x0.1 degrees on long-lat projection,
  and 34 vertical levels.

Download the input via the catalog tool (:numref:`sec-ModelCode`) as follows:

.. code-block:: bash

    # download 2018 meteorology for the EMEP01 domain
    catalog.py -Y 2018 -m --met-domain EMEP0201

    # download other input files
    catalog.py --input

The meteorology files will be placed under
``EMEP_MSC-W_model.v5.5.OpenSource/meteo2018/EMEP0201/``,
and the remaining input files will be placed under
``EMEP_MSC-W_model.v5.5.OpenSource/input/``

This are all input files needed to run the EMEP/MSC-W model,
except the aircraft emissions (``AircraftEmis_FL.nc``),
and forest fire emissions (``FINN_ForestFireEmis_2018.nc``).
See sections :numref:`emisair` and :numref:`emisff`
for details about these emissions data.

IMPORTANT:
    The input data available in the EMEP/MSC-W Open Source Web site should
    be appropriately acknowledged when used for model runs. If nothing else
    is specified according to references further in this chapter, please
    acknowledge EMEP/MSC-W in any use of these data.


.. csv-table:: List of input data files
    :name: tab-inputdata
    :header: **Data**, **Name**, **Format**
    :delim: &

    **Meteorology data**& ``meteoYYYY/GRID``&
    Meteorology       & ``meteoYYYYMMDD.nc`` (365+1 files)            & netCDF [#YMD]_
    Degree-day factor & ``DegreeDayFactors.nc``                       & netCDF
    **Other Input files**& ``input/``&
    Global Ozone      & ``Logan_P.nc``                                & netCDF [#O3]_
    BVOC emissions    & ``EMEP_EuroBVOC.nc``                          & netCDF
    Landuse           & ``glc2000xCLMf18.nc`` and ``Landuse_PS_5km_LC.nc``& netCDF
    N depositions     & ``AnnualNdep_PS50x_EECCA2005_2009.nc`` or ``CAMS-GLOB-SOIL_Glb_0.5x0.5_soil_nox_v2.4clim_monthly.nc`` & netCDF
    Road dust         & ``RoadMap.nc`` and ``AVG_SMI_2005_2010.nc``   & netCDF [#Optional]_
    Aircraft emissions& ``AircraftEmis_FL.nc``                        & netCDF [#Optional]_
    Surface Pressure  & ``SurfacePressure.nc``                        & netCDF [#Optional]_
    Forest Fire       & ``FINN_ForestFireEmis_YYYY.nc``               & netCDF [#Optional]_
    Dust files        & ``Soil_Tegen.nc``                             & netCDF [#Optional]_
                      & ``SoilTypes_IFS.nc``                          & netCDF [#Optional]_
    Emissions         & ``EECCA/emislist.POLL`` (7 files, EMEP 50km PS grid)             & ASCII [#POLL]_
                      & ``EMEP01/GNFRemis_EMEP01_2015.nc`` (regional, :math:`0.1\times 0.1`  lon-lat) & netCDF [#Optional]_
    Vertical level distribution         & ``Vertical_levels20_EC.txt`` (for 20lev runs) & ASCII
    Time factors for monthly emissions  & ``MonthlyFac.POLL`` (7 files)       & ASCII [#POLL]_
    Time factors for daily emissions    & ``DailyFac.POLL`` (7 files)         & ASCII [#POLL]_
    Time factors for hourly emissions   & ``HourlyFacs.INERIS``               & ASCII
    Natural |SO2|                       & ``DMS.nc``                          & netCDF
    Volcanoes                           & ``columnsource_emission.csv``       & ASCII
                                        & ``columnsource_location.csv``       & ASCII
    Lightning emissions                 & ``lt21-nox.datMM`` (12 files)       & ASCII [#YMD]_
    Emissions speciation                & ``emissplit.defaults.POLL``         & ASCII [#POLL]_
                                        & ``emissplit.specials.POLL``         & ASCII [#POLL]_ [#Optional]_
    Emission factors for scenario runs  & ``femis.dat``                       & ASCII
    Photo-dissociation rates            & ``jclear.SEASON`` (4 files)         & ASCII [#SEASON]_
                                        & ``jcl1.SEASON`` (4 files)           & ASCII [#SEASON]_
                                        & ``jcl3.SEASON`` (4 files)           & ASCII [#SEASON]_
    Landuse definitions                 & ``Inputs_LandDefs.csv``             & ASCII
    Stomatal conductance                & ``Inputs_DO3SE.csv``                & ASCII
    Sites locations for surface output  & ``sites.dat``                       & ASCII
    Sondes locations for vertical output& ``sondes.dat``                      & ASCII

.. rubric:: Footnotes
.. [#YMD] ``YYYY``: year, ``MM``: month, ``DD``: day.
.. [#O3] |O3| boundary condition data in 30 levels.
.. [#Optional] Optional, in most cases.
.. [#POLL] ``POLL``: pollutant type (|NH3|\ , CO, |NOx|\ , |SOx|\ , NMVOC, |PM25| and |PMco|\ ).
.. [#SEASON] ``SEASON``: seasonal files (jan, apr, jul, oct).


NetCDF files
------------

Meteorology
~~~~~~~~~~~

The daily meteorological input data (``meteoYYYYMMDD.nc``, where ``YYYY`` is
year, ``MM`` is month and ``DD`` is day) used for the EMEP/MSC-W Model are based
on forecast experiment runs with the Integrated Forecast System (IFS), a
global operational forecasting model from the European Centre for
Medium-Range Weather Forecasts (ECMWF).

The IFS forecasts has been run by MSC-W as independent experiments on
the HPCs at ECMWF with special requests on some output parameters.
The meteorological fields are retrieved on a
:math:`0.1^\circ\times 0.1^\circ` longitude latitude coordinates and
interpolated to :math:`0.3^\circ\times 0.2^\circ`.
Vertically, the fields on 60 eta (\ :math:`\eta`\ ) levels from the IFS model are
interpolated onto the 37 EMEP eta (\ :math:`\eta`\ ) levels. The meteorology is prepared
into 37 eta levels since the model is under test for a finer vertical resolution.

The open source code is released with 20 eta levels and
to make the model read the meteorology properly, a description of the 20
vertical sigma levels is needed. This is provided in an ASCII file
called ``Vertical_levels20_EC.txt`` together with the other input data (:numref:`tab-inputdata`).
The version of the IFS model used for preparing these fields for 2018 and earlier years,
Cycle 40r1, is documented in
https://www.ecmwf.int/en/forecasts/documentation-and-support/changes-ecmwf-model/cycle-40r1/cycle-40r1.
2019 and later years are based on Cycle 48r1, described in
https://confluence.ecmwf.int/display/FCST/Implementation+of+IFS+Cycle+48r1.
Some verification and description of 2018
meteorological fields are given in Chapter 2 of the EMEP Status Report 1/2020
https://www.emep.int/mscw/mscw_publications.html#2020.

Acknowledgement:
    ECMWF, met.no


.. csv-table:: Input meteorological data used in the EMEP/MSC-W Model
    :name: tab-metinput
    :header: **Parameter**, **Unit**, **Description**
    :delim: &

    **3D fields** && for 37 :math:`\eta`
    :math:`u_wind, v_wind`              & :math:`m/s`       & Horizontal wind velocity components
    :math:`etadot``                     & :math:`Pa/s``     & Vertical velocity in :math:`\eta` coords
    :math:`specific_humidity`           & :math:`kg/kg`     & Specific humidity
    :math:`potential_temperature`       & :math:`K`         & Potential temperature
    :math:`cloudwater`                  & :math:`kg/kg`     & Cloud water
    :math:`cloudice`                    & :math:`kg/kg`     & Ice cloud water
    :math:`3D_cloudcover`               & :math:`\%`        & 3D Cloud cover
    :math:`convective_updraft_flux`     & :math:`kg/m^2/s`  & Convective updraft flux
    :math:`convective_downdraft_flux`   & :math:`kg/m^2/s`  & Convective downdraft flux
    :math:`precipitation`               & :math:`kg/m^2`    & Precipitation
    **2D fields** && for surface
    :math:`surface_pressure`            & :math:`hPa`       & Surface pressure
    :math:`temperature_2m`              & :math:`K`         & Temperature at :math:`2 m` height
    :math:`relative_humidity_2m`        & :math:`\%`        & Relative humidity at :math:`2 m` height
    :math:`surface_flux_sensible_heat`  & :math:`W/m^2`     & Surface flux of sensible heat
    :math:`surface_flux_latent_heat`    & :math:`W/m^2`     & Surface flux of latent heat
    :math:`surface_stress`              & :math:`N/m^2`     & Surface stress
    :math:`sea_surface_temperature`     & :math:`K`         & Sea surface temperature
    :math:`soil_water_content`          & :math:`m^3/m^3`   & Soil water content
    :math:`deep_soil_water_content`     & :math:`m^3/m^3`   & Deep soil water content
    :math:`large_scale_precipitations`  & :math:`m/s`       & Large scale precipitation
    :math:`convective_precipitations`   & :math:`m/s`       & Convective precipitation
    :math:`snow_depth`                  & :math:`m`         & Snow depth
    :math:`fraction_of_ice`             & :math:`\%`        & Fraction of ice
    :math:`SMI1`                        &                   & Soil moisture index level 1
    :math:`SMI3`                        &                   & Soil moisture index level 3
    :math:`pblh``                       & :math:`m`         & Planetary boundary layer height
    :math:`u10, v10`                    & :math:`m/s`       & Wind at :math:`10 m` height

.. _`emisnew`:

Gridded emissions
~~~~~~~~~~~~~~~~~

Gridded emissions in NetCDF care be used in conjunction with sector definitions.
See section ``Defining own sectors`` and ``Country Variable (CV) format``

    The main advantage of NetCDF emissions is that all the information
    about the data (projection, units) is given in the same file. This
    allows the code to reproject the emissions to any grid projection on
    the fly. It is easy to visualize the emissions of one country with
    simple tools, like ncview. The data is simple to interpret and it is
    possible to add new countries to an existing file (with appropriate
    tools).


Global Ozone
~~~~~~~~~~~~

Initial concentration of ozone are required in order to initialize the
model runs. Boundary conditions along the sides of the model domain and
at the top of the domain are then required as the model is running.

The ``Logan_P.nc`` file contains monthly averaged fields in NetCDF format.
The initial and background concentrations are based on the Logan (1998)
climatology. The Logan climatology is scaled on run time according to the
Mace Head measurements as described in Simpson *et al.* (2003). For a
number of other species, background/initial conditions are set within
the model using functions based on observations (Simpson *et al.*, 2003
and Fagerli *et al.*, 2004).

BVOC emissions
~~~~~~~~~~~~~~

Biogenic emissions of isoprene and monoterpene are calculated in the
model as a function of temperature and solar radiation, using the
landuse datasets. The light and temperature dependencies follow
the ideas proposed
in Guenther et al (1993,1995), the first step in the
emission processing is to define 'standard' emission potentials, which
give the emissions of particular land-covers at standard environmental
conditions (:math:`30^\circ C`  and photosynthetically active radiation of
1000 :math:`\mu mole/m^2/s`).

European forests are treated in most detail.  For these, BVOC emission
potentials have been created from the the map of forest species generated
by Koeble and Seufert (2001). This work provided maps for 115 tree
species in 30 European countries, based upon a compilation of data from
the ICP-forest network. The emission potentials for each species are
as given in Simpson *et al.*, 2012, and have been aggregated into the
four default forest classes used by EMEP over Europe (DF, CF, NF, BF).
The NetCDF file ``EMEP_EuroBVOC.nc`` provides the aggregated emission
potenitals for these 4 categories. These emission potentials have unit
:math:`\mu g/m^2/h`\ , and refer to emissions per area of the appropriate
forest category.

On the global scale, new landcover maps were created as a combination of
GLC2000 and Community Land Model (CLM) data as described in 
Simpson *et al.*, 2017.
The default emission potentials
are given for these extra CLM categories, and for any non-forest land-cover on
Europe  in the file
``Inputs_LandDefs.csv``. The underlying emission potentials, land-cover
data bases, and model coding have however changed substantially since
model version v.2011-06. The new approach is documented in Simpson *et
al.*, 2012 and Simpson *et al.* 2017.

Landuse
~~~~~~~

Landuse data are required for modelling boundary layer processes (i.e.
dry deposition, turbulent diffusion). The EMEP/MSC-W model can accept
landuse data from any data set covering the whole of the domain,
providing reasonable resolution of the vegetation categories. Gridded
data sets providing these landuse categories across the EMEP domain have
been created based on the data from the Stockholm Environment Institute
at York (SEI-Y) and from the Coordinating Center for Effects (CCE). 16
basic landuse classes have been identified for use in the deposition
module in the model, and three additional "fake" landuse classes are
used for providing results for integrated assessment modeling and
effects work.

There are two NetCDF files included, one file
``Landuse_PS_5km_LC.nc`` on 5 km resolution over the EMEP domain,
and a global ``LanduseGLC.nc`` which combines data from GLC2000 with the Community Land Model (CLM).
The different landuse types are desribed
in Simpson et al (2012) and Simpson et al. (2017).

Degree-day factor
~~~~~~~~~~~~~~~~~

Domestic combustion which contribute to a large part of GNFR C (SNAP 2), varies
on the daily mean temperature. The variation is based on the heating
degree-day concept. These degree days are pre-calculated for each day
and stored in the file ``DegreeDayFactors.nc``. See Simpson et al. (2012)
section 6.1.2.

|NOx| depositions
~~~~~~~~~~~~~~~~~

Areas with high NO deposition loads have greater soil-NO emissions. To
include this in the model, a NetCDF file where pre-calculated
N-depositions are included. The file made by the results from the
EMEP/MSC-W model runs over a 5-year period.

Road Dust
~~~~~~~~~

Road traffic produces dust. These can be handled as separate emission inputs in
the EMEP/MSC-W model using the ``Emissions_mod.f90`` module. To include road
dust, set ``USE_ROADDUST=.true.`` in ``config_emep.nml``. There are two
files included in input data, ``RoadMap.nc`` and ``AVG_SMI_2005_2010.nc``.
``RoadMap.nc`` include gridded roads and PM emissions over Europe, while the
Soil Moisture Index (SMI) file ``AVG_SMI_2005_2010.nc`` used to estimate
emissions is global. Hence road dust emissions can currently only be calculated
for the European domain. However, some countries for which road dust is important
(e.g., Scandinavian countries), reported emissions already include road dust. By
default we therefore set ``USE_ROADDUST=.false.``, with road dust as a separate
emission source effectively being deprecated. 

.. _`emisair`:

Aircraft emissions
~~~~~~~~~~~~~~~~~~

In the EMEP/MSC-W model aircraft emissions are 'ON' by default.
It can be switched 'OFF' by setting ``USES%AIRCRAFT_EMIS=F`` in ``config_emep.nml``.
When using aircraft emissions there are two options:
Either the dataset provided by the EU-Framework Programme 6 Integrated
Project QUANTIFY (http://www.pa.op.dlr.de/quantify), or the more recent
CAMS-GLOB-AIR dataset which can be downloaded from the Atmosphere Data Store
(https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-global-emission-inventories?tab=form).
However, before using these data a protocol has to be signed, which is why
the data file can not be provided directly on the EMEP MSC-W Open Source
website.
Since rv4.39, CAMS-GLOB-AIR has been the default dataset. Registering
at the Atmosphere Data Store is straightforward
(https://ads.atmosphere.copernicus.eu/user/register?destination=/cdsapp)
(If you rather want to use older QUANTIFY dataset go to
http://www.pa.op.dlr.de/quantify, click on
'QUANTIFY emission inventories and scenarios', and then on
'Register'. That page will provide information about the registration
process and the protocol that has to be signed. Once you are registered,
click 'Login' and provide user name and password. On the new page,
search for 'Emissions for EMEP', which links directly to the ``Readme`` file
and the emission data file in NetCDF format. Download the emission data
file and place it in the input folder.)

Natural |SO2|
~~~~~~~~~~~~~~~~~~~~

Natural |SO2| emissions (dimethylsulfide (DMS) from sea) are
provided as monthly gridded files. The values are computed taking into account sea surface
temperature and wind speed. Surface water concentrations of DMS (needed for the flux
calculation) are taken from SOLAS (Surface Ocean Lower Atmosphere Study) and were
downloaded from https://www.bodc.ac.uk/solas_integration/implementation_products/group1/dms/.

Surface Pressure
~~~~~~~~~~~~~~~~

If ``USE_AIRCRAFT_EMIS=.true``. in ``config_emep.nml``, then in
addition to the Aircraft Emission file, there will be need for a
``SurfacePressure.nc`` file, which is already in the ``/input`` folder. The
NetCDF file consists of surface pressure fields for each of the months
in 2008 called ``surface_pressure``, and one field for the whole year
called ``surface_pressure_year``. All fields are given in Pa.

.. _`emisff`:

Forest Fire
~~~~~~~~~~~

Since model version rv3.9 (November 2011), daily emissions from forest
and vegetation fires are taken from the "Fire INventory from NCAR
version 1.0" (FINNv1, Wiedinmyer et al. 2011). Data are available from
2005, with daily resolution, on a fine :math:`1 km\times1 km` grid.
We store these data on a slightly coarser grid (\ :math:`0.2^\circ\times 0.2^\circ`\ )
globally for access by the EMEP MSC-W model. To include forest fire
emissions set ``USE_FOREST_FIRES=.true.`` in ``config_emep.nml`` and
download the 2012 GEOS-chem daily data
http://bai.acom.ucar.edu/Data/fire/. The data needs to be stored with
units mole/day in a NetCDF file called ``FINN_ForestFireEmis_2015.nc``
compatible with the ``ForestFire_mod.f90`` module.

Dust files
~~~~~~~~~~

The annual ASCII data for sand and clay fractions as well as the monthly
data for boundary and initial conditions for dust from Sahara are
replaced with a single NetCDF file ``Soil_Tegen.nc`` since 2013. This
covers data for a global domain in :math:`0.5\times 0.5` degree
resolution.

The variables 'sand' and 'clay' gives the fraction (in %) of sand an
clay in the soil for each grid cell over land.

The files are used by the module ``DustProd_mod.f90``, which calculates
windblown dust emissions from soil erosion. Note that the
parametrization is still in the development and testing phase, and is by
default 'turned off'. To include it in the model calculations, set
``USE_DUST=.true.`` in ``config_emep.nml``. The user is recommended to
read carefully documentation and comments in the module ``DustProd_mod.f90``.

There is also a possibility to include boundary and initial conditions
for dust from Sahara. The input file gives monthly dust mixing ratios
(MM - month, e.g. 01, 02, 03,...) for fine and coarse dust from Sahara.
The files are based on calculations from a global CTM at the University
of Oslo for 2000. To include Saharan dust, set ``USE_SAHARA=.true.`` in
``config_emep.nml``.

Another source for dust is an arid surface. This is accounted for by
soilmosture calculations in ``DustProd_mod.f90``. Together with Soil
Water Index from the meteorology files and permanent wilting point (pwp)
from ``SoilTypes_IFS.nc``. This file is global and NetCDF. See Simpson et
al. (2012) section 6.10.

ASCII files
-----------


Volcanoes
~~~~~~~~~

Emissions from volcanic passive degassing of |SO2| are included
for the active Italian volcanoes, Etna, Vulcano and Stromboli, and based upon the
officially submitted data. To consider these volcanic emissions, we need
to feed the locations and heights of volcanoes into the model. The input
file ``columnsource_location.csv`` contains the geographical coordinates
(latitudes and longitudes) and the heights (in meters) of the included
volcanoes, while ``columnsource_emission.csv`` contains the emission
parameters.

Since 2010 the EMEP/MSC-W  model has also been used to model the transport of
ash and |SO2| from volcanic eruptions. In addition to data for
passive degassing of |SO2|\ , the above two input files also
contain locations and emission parameters for two recent eruptions of
Icelandic volcanoes (Eyjafjallajökull in 2010 and Grimsvötn in 2011).
In order to include emissions from these eruptions one needs to set
``USE_ASH=.true.`` in ``config_emep.nml``.

Gridded emissions
~~~~~~~~~~~~~~~~~

The official emission input for the EMEP/MSC-W model consists of gridded
annual national emissions based on emission data reported every year to
EMEP/MSC-W (until 2005) and to CEIP (from 2006) by each participating
country. More details about the emission input with references can be
found in Chapter 4 of the EMEP Status Report 1/2003 Part I (Simpson et
al., 2003).

Since 2015 different formats of gridded emissions can be used and mixed
(with some restrictions) in the EMEP model under one common framework.
The new emission system is described in :numref:`emisnew`. Here we focus
only on the "old style" ASCII emission format.

Seven gridded emission input files (``emislist.poll``) are available in
ASCII format for the following compounds: CO, |NH3|\ ,
|NOx|\ , |PM25|\ , |PMco|\ , |SOx| and VOC.

The gridded ASCII emission files contain 16 columns where the first
column represents the country code
(http://www.emep.int/grid/country_numbers.txt), the second and the third
columns are the :math:`i` and :math:`j` indices of the EMEP grid, the fourth and
fifth columns are the total emissions from low and high sources, and the
last 11 columns contain emissions from 10 anthropogenic SNAP sectors
(http://reports.eea.eu.int/technical_report_2001_3/en) and 1
source-sector called "Other sources and sinks", which include natural and
biogenic emission sources. The data are given with the :math:`Mg`\ .

Acknowledgement:
    EMEP

Time factors for emissions
~~~~~~~~~~~~~~~~~~~~~~~~~~

Monthly and daily time factors for emission can be specified for the 7
compounds (CO, |NH3|\ , |NOx|\ , |PM25|\ , |PMco|\ , |SOx| and VOC).  There is
one file available per compound in ASCII format.

The first two columns in the files represent the country code
(http://www.emep.int/grid/country_numbers.txt), the second column represents an
index that can be referenced by the sector definion (orginally this index
corresponded to a SNAP sector). 

In the monthly files (``MonthlyFacFile``), the 12 consecutive columns represent
the time factors corresponding to the months of the year. These time factors
are interpolated according to whether the "current" simulation day is in the
first or second half of the month. For days in the first half of the month, the
monthly factor is a combination of the factor of the previous month and the
current month. For days in the second half of the month, the monthly factor
combines the factor of the current and the next month (for details, see the
source code).

In the daily files (``DailyFacFile``) there are 7 consecutive columns
representing the time factor for each day of the week. The monthly timefactors
are normalized in such a way that when combined with the daily timefactors, the
total emission stays the same.

The file defined in ``HourlyFacFile`` includes factors for each of the eleven
SNAP sectors for every hour (the columns) for each day of the week, see Simpson
et al. (2012) section 6.1.2. An additional file defined in
``HourlyFacSpecialsFile`` can be created by the user with modified hourly
factors to be used for specific countries. The format is the same as for the
default factors, except for an additional first column speicifying the country
code number.


Emission heights
~~~~~~~~~~~~~~~~

In previous versions the emission height distributions was given in a separate file. Now it is 
part of the code, and can also be modified by the users using config_emep.nml setting
(see section "defining own sectors).

A set of vertical distribution for different sectors are predefined in the model. 
The release heights are defined as a set of fractions released into predefined layers.

The release height definitions are independent of the layers used by the model. 

There are 8 predefined release heights distributions. Those can also be defined through the config_emep.nml setting. The following will give exactly the same distributions as the predefined. You can then modify the values, or add new defined distributions.

.. code-block:: Fortran
    :caption: Default definition of emission height distributions

    Emis_Zlevels(1:)20.0,   50.0,   92.0,  184.0,  324.0,  522.0,  781.0, 1106.0,
    Emis_h(1:,1) = 0.000,  0.000,  0.000,  0.003,  0.147,  0.400,  0.300,  0.150, 
    Emis_h(1:,2) = 1.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
    Emis_h(1:,3) = 0.060,  0.067,  0.093,  0.750,  0.030,  0.000,  0.000,  0.000, 
    Emis_h(1:,4) = 0.050,  0.063,  0.087,  0.700,  0.100,  0.000,  0.000,  0.000, 
    Emis_h(1:,5) = 0.020,  0.034,  0.046,  0.600,  0.300,  0.000,  0.000,  0.000,
    Emis_h(1:,6) = 0.000,  0.000,  0.000,  0.410,  0.570,  0.020,  0.000,  0.000,
    Emis_h(1:,7) = 0.200,  0.300,  0.020,  0.044,  0.066,  0.094,  0.123,  0.153, 
    Emis_h(1:,8) = 0.200,  0.800,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,

The ``Emis_Zlevels`` defines the height of the layer boundaries for emissions in meters. (Standard atmosphere is assumed to transform those in Pressure by the model). The first layers is from surface to 20 meters, the second layer from 20 to 50 m... until the eigth and last layer which runs from 781 to 1106 meters. 

For example sectors defined with the height index "1", will release nothing in the three lowest layers, 0.3% into the fourth layer, 14.7% into the fifth layer etc. 

The layers defined in Emis_h are independent from the layers used in the model run and do not need to be adapted if the number of model layers is modified. 
The actual resulting distribution of emissions into model layers is computed by the model and will be shown in the standard output.



.. _`sec-femis`:

Emission factor for scenario runs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Scenario run in the case of the EMEP/MSC-W model means a run
to test the impact of one or more pollutants from a particular country.

Emission factors are applied to specified countries and emission sectors
and can be set by changing the ASCII file ``femis.dat``. This file can
be changed by the users according to their needs.
See the `Source Receptor (SR) Runs <https://emep-ctm.readthedocs.io/en/latest/Subrun.html#source-receptor-sr-runs>`_ section for details.


Vertical_levels.txt
~~~~~~~~~~~~~~~~~~~

Defines the vertical model layers. The numbers in Vertical_levels.txt correspond to the "A" and "B" coefficients of the hybrid (eta) coordinates (P=A+B*Psurf).

Close to the surface, A should be small, and higher up we should use pressure levels. Then there is a gradual transition from surface to pressure levels.

If the file is not provided, the meteorological vertical levels are used. In principle the model levels can be completely different, but it is more sensible to define layer boundaries that match the meteorological levels.

NB: it is important not to define a lowest layer thinner than about 45 meters; the deposition scheme will fail if the middle of the lowest layer is smaller than the highest defined vegetation.

Chemical speciation of emissions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Many of the emission files give emissions of a group of compounds, e.g.
|NOx| includes NO+|NO2|\ , and VOC can include many compounds. The
information needed to retrieve emissions of individual compounds from
these the gridded files is given in files labelled
``emissplit.defaults.POLL`` or ``emissplit.specials.POLL``,
where ``POLL`` can be |NOx|\ , VOC, etc.

The defaults file give the emission split for each sector split index 
(one per row, with second index being the sector split index), which is applied to all
countries by default. For VOC this split was derived from the UK
inventory of Passant (2002), as part of the chemical comparison project
of Hayman *et al.* (2011).

The specials files are in general optional, and can be used to specify
speciation for particular countries or sectors. The 1\ :sup:`st`
column specifies the country code of interest, the second the
sector index.

If forest fires are used, then the file ``emissplit.specials.voc`` is
required (not optional), and the country-code 101 used to specify the
VOC speciation of forest fires in this file.

Lightning emissions
~~~~~~~~~~~~~~~~~~~

Emissions of |NOx| from lightning are included in the model as
monthly averages on T21 (\ :math:`5.65^\circ\times 5.65^\circ`\ )
resolution (Køhler *et al.*, 1995). The lightning emissions are defined
on a :math:`64\times 32` grid with 17 vertical levels, with global
coverage, and are provided as 12 ASCII files ``lightningMM.dat``.

Landuse definitions
~~~~~~~~~~~~~~~~~~~

For the vegetative landuse categories where stomatal modelling is
undertaken, the start and end of the growing season (SGS, EGS, in days)
must be specified. The calculation of SGS and EGS with respect to
latitude is done in the module ``LandDefs_mod.f90``. The parameters
needed to specify the development of the leaf area index (LAI) within
the growing season are given in the ASCII file ``Inputs_LandDefs.csv``.
For more information, see chapter 5 of the EMEP Status Report 1/2003
Part I (Simpson *et al.*, 2003).

The file, designed to be opened with excel or gnumeric, contains a
header briefly explaining the contents of the 14 columns. The first
three columns are representing the landuse name, code (which are
consistent with those in ``Landuse.Input`` file) and type (grouping of the
landuse classes). The fourth column (PFT) gives a plant-functional type
code (for future use), the fifth gives the maximum height of vegetation
(m), the sixth indicates albedo (%) and the seventh indicates
possible source of |NHx| (0 off/1 on, curently not used).
Columns 8 to 11 define the growing season (day number), column 12 and 13
lists the LAI minimum and maximum (\ :math:`m^2/m^2`\ ) and columns 14
and 15 defines the length of the LAI increase and decline periods (no.
of days). Finally, the last four columns give default values of foliar
biomass and biogenic VOC emission potentials. See Simpson et al., (2012)
for details.

Stomatal conductance
~~~~~~~~~~~~~~~~~~~~

Parameters for the stomatal conductance model, deposition of
|O3| and stomatal exchange (DO3SE) must be specified. That are
based upon the ideas in Emberson *et al.*, 2000, and are discussed in
Simpson and Emberson, 2006 and Tuovinen et al. 2004.

The ASCII file ``Inputs_DO3SE.csv`` provides land-phenology data of each
landuse type for stomatal conductance calculations. The data are
summarised in Table 5.1 in Chapter 5 of the EMEP Status Report 1/2003
Part I (Simpson *et al.*, 2003).

The file contains a **header** with the contents of the file, with
different factors needed for each of the landuse classes used in the
EMEP/MSC-W model. The first two columns represent the landuse code
(which are consistent with those in ``Landuse.Input`` file) and name.
The next 22 values are different phenology factors.

Photo-dissociation rates
~~~~~~~~~~~~~~~~~~~~~~~~

The photo-dissociation rates (J-values) are provided as lookup tables.
The method is previously described in Jonson *et al.*, (2001). J-values
are provided as clear sky, light cloud and dense cloud conditions, and
the model interpolates between these according to cloudiness from the
meteorological input data. In the lookup tables data are listed for
every 10 degree latitude at an interval of 1 degree zenith angle at
every model height.

For the two types of cloud conditions there are one ASCII file
averaged for each season (``SS``); 01, 02, 03 and 04. For light cloud the
four seasonal files are called ``jcl1kmSS.dat``, for dense cloud
conditions the four seasonal files are called ``jcl3kmSS.dat``, and then
for clear sky four files called ``jclearSS.dat``. In addittion
there are two files for June called ``jcl1.jun`` and ``jcl3.jun``.

Each file contains 18 columns. The first column is latitude of zenith
angle and then the next 17 are the values for the model levels with the
1/s. For more details about these rates, please read
Chapter 7.2 of the EMEP Status Report 1/2003 Part I (Simpson *et al.*,
2003).

.. _`sec-sitessondes-input`:

Site and Sonde files
~~~~~~~~~~~~~~~~~~~~

The model provides a possibility for
extra output data of surface concentration for a set of specified
measurement site locations and concentrations for the vertical column
above a set of specified locations. These site and sonde locations are
listed in the ASCII files ``sites.dat`` and ``sondes.dat``
files. These files can be changed by the user, and provide the outputs
described in :numref:`sec-sitesonde`.

Two main options are available for the output of ASCII files for
comparison with measurements or detailed model analysis. These are

sites
    output of surface concentrations for a set of specified measurement
    site locations.

sondes
    output of concentrations for the vertical column above a set of
    specified locations.

Both sites and sondes are specified and handled in similar ways, in the
module ``Sites_mod.f90``, so we treat them both together below.
Locations are specified in the input files ``sites.dat`` and ``sondes.dat``.
The files start with a description of its content followed by a list of
the stations. For example, a sondes.dat input file may look like this:
.. literalinclude:: sites.dat
    :caption: Site location definition (``sites.dat``) example.

.. literalinclude:: sitesLLHrel.dat
    :caption: Site location definition, LatLonHrel Coords example.

.. literalinclude:: sitesLLZ.dat
    :caption: Site location definition, LatLonZm Coords example.

The first line in each file is a header with file content. Then, the
contents are described in more detail. Text strings after ``#`` are just
clarifying comments. 'Area', e.g., is the domain to which the stations
belong, e.g. 'Northern Hemisphere'.

Text after ``:`` is read in by the model:

Units
    Either 'deg' (degrees) or 'index' (model grid indices).

Coords
    Vertical coordinate system that is used in the model - see below.

The VertCoords system was changed in rv4.36, to allow several options.
Specifying altitude rather than model coordinate is of course an obvious
alternative to the model layer number that was previously used, but it is also easily mis-used and mis-understood too. For example, a site at 500m might need vertical profile and/or deposition correction if sitting on a 500m high plateau (hence it has a relative altitude of 0m, and should use iz=KMAX_MID), but if sitting on an isolated mountain in terrain of height 0m, the relative altitude would indeed to 500m and we should pick from some model level which represents that. Tricky! (Best might be to give station altitude and relative altitude in separate columns, so it is explicit at least.)

Briefly though, one can now set the "Coords" parameter (in the header of the sites.dat input file, used to be just "LatLong") to be one of:

  1. LatLonKdown  - same as older LatLong, with lat/lon coordinates, then k-number with ground level being e.g. 20

  2. LatLonZm - give lat, long, then altitude in metres. The model will then compare that altitude with the model's topography, to estimate a relative altitude (Hrel). It then calculates which model layer corresponds to that altitude.

  3. LatLonHrel - give lat, long, and then your own preferred relative altitude (Hrel). This relative altitude may be calculated by comparison to digital elevation model (see for example, Loibl, W.; Winiwarter, W.; Kopsca, A.; Zufger, J. & Baumann, R. Estimating the spatial distribution of ozone concentrations in complex terrain Atmos. Environ., 1994, 28, 2557-2566).

  4. IJKdown - does everything in model's i, j and k coordinates

In principle (3) should produce the best results, with relative altitudes calculated compared to local topography (e.g. < 5km). Option (3) also allows the same Hrel to be used regardless of the underlying meteo resolution. One can also get Hrel from the TOAR database for some sites (Schultz, MG, et al 2017 Tropospheric Ozone Assessment Report: Database
and metrics data of global surface ozone observations. Elem Sci Anth,
5: 58, DOI: https://doi.org/10.1525/elementa.244).

However, testing of these systems against e.g. diunral profiles of ozone shows rather unpredictable results in some cases, and often just using the model's surface concentrations procudes results which are as good as those based upon altitude. The new systems are very fliexble thouhg, and allow the user to explore different methodologies.

Both ``sites.dat`` and ``sondes.dat`` files are optional, but recommended.
The species and meteorological data requested for site and sonde output are
specified in ``Config_module.f90`` by the use of arrays.
Only a few met fields are defined so far but more can be added into
``Sites_mod.f90`` as required.

The output files ``sites_2015.csv`` and ``sondes_2015.csv`` are comma
separated files that can be read by excel. netcdf versions of these files
are also provided, ``sites_2015.nc`` and ``sondes_2015.nc``.
If you include the whole year, or the 31\ :sup:`st` December,
``sites_2016.csv`` and ``sondes_2016.csv`` are also included in the output.


