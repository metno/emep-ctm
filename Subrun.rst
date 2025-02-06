.. _`ch-submitarun`:

Setting the input parameters
============================

In this chapter we provide detailed information on how to configure the parameters of the regional EMEP/MSC-W model.

In general the parameters and pathes to the input files can all be set in the configuration file ``config_emep.nml`` (a fortran namelist).

``config_emep.nml``
-------------------

The default parameter, constants and flags are defined in ``Config_module.f90``, and they can be overwritten by ``config_emep.nml`` settings. 

Here is an example of content of the most important parameters:

.. code-block:: text
  :name: config-emep
  :caption: Basic namelist example; ``config_emep.nml`` extract.

  &Model_config
    GRID = 'EECCA',
    iyr_trend = 2015,
    runlabel1 = 'Base',
    runlabel2 = 'Opensource_Setup_2022',
    startdate = 2015,01,01,00,
    enddate = 2015,01,10,24,
    spinup_enddate = 2015,01,02,06,

    DataPath(1) = '../input', ! define 'DataDir' keyword

    meteo = '../meteoYYYY/GRID/meteoYYYYMMDD.nc',
    DegreeDayFactorsFile = 'MetDir/DegreeDayFactors.nc',
  !------------------------------
    EmisDir = 'DataDir/GRID',
    emis_inputlist(1)%name= 'EmisDir/gridPOLL', !example of ASCII type
    emis_inputlist(1)%type = 'GNFR_CAMSsectors',
  !--------Sub domain x0, x1, y0, y1
    RUNDOMAIN = 36, 100, 50, 150, ! EECCA sub-domain
  &end

In the extract above, the model is run for the period 1 January to 10 January 2015 and the trend year used is 2015,
with output staring at 6 UTC on January 2nd.
The ``startdate`` must correspond to a time step on the input meteorology,
``enddate`` and ``spinup_enddate`` can be chosen freely.
Output files will be stored with the name 'Base' and the meteorological correspond to the 'EECCA' grid.

Note the definition of the meteorological input ("meteo"):
the keywords, 'YYYY', 'MM', 'DD' and 'GRID' will be replaced on the fly by respectively, year, month, day and grid ('EECCA' here).

It is possible to run the model on a smaller domain than the full
regional model domain, as defined by indexes :math:`x` and :math:`y`.
For the 'EECCA' grid  :math:`x=1,\ldots,132; y=1,\ldots,159`\ .
To set a smaller domain, use ``RUNDOMAIN`` variable in the ``Model_config``
namelist to indicate the sub-domain indexes. In the config_emep extract above,
``RUNDOMAIN`` defines a subdomain with :math:`x=36,\ldots,100; y=50,\ldots,150`\ .
The indices refer always to the meteo grid (starting from 1).
If the indices are outside the grid, only the portion within the full grid is taken into account. 


Base run
--------

To run the model you need access to the executable (``emepctm``) and the ``config_emep.nml`` files.
The path to the other input files are defined by defaults or set in the ``config_emep.nml`` file.

To run the model you can either run interactively (mostly for short run), using something similar to (it depends on your system):

``mpirun emepctm``

For longer runs you should run the model as a batch job (the details will depend on your system).
If the run was successful, a message

.. code-block:: text

    programme is finished

will be stated at the end of the standard output (screen or some appropriate file for batch jobs).

The model results will be written to the same directory. Please make sure there is enough disk place for the model
results. For more information about the model result output files,
see :numref:`ch-output`.

If for some reason the model crashed, please check both the log and the
error file for any clue of the crash. After fixing the problem the job
can be submitted again. 

The variables wanted in the output are specified in the
``OutputConcs``, ``DDEP_ECOS``, ``DDEP_WANTED``, ``WDEP_WANTED`` and ``OutputVegO3``
parameters for surface concentrations, depositions and other miscellaneous outputs.

The main ouput files are
 - Base_fullrun.nc gives values averaged over the entire simulation.
   This file has always 1 record and the data is an average over the the period startdate to enddate. 
 - Base_month.nc average over each calendar month (12 records for a yearly run) 
 - Base_day.nc averages over 24 hours, starting at 06:00 UTC (each record gives average values between 06:00 and 06:00 the next day.
   The first time step is written out at 06:00 hr on 2Jan (i.e the first 6 hours are “lost”).
   The last record may also be averaged only until the end of the run, and not at 06:00.)
 - Base_hour.nc, averaged over 1 hour.
 - Base_hourInst.nc, instantaneous values every hour.


Source Receptor (SR) Runs
-------------------------

The EMEP/MSC-W model can be used to test the impact of reduced emission
of one or more pollutants from a particular country or a number of
countries. Such runs are called "Scenario runs". They are used for source-receptor calculations.

Emission factors for reduced emissions of pollutants from different
sectors and countries can be defined in the input file ``femis.dat``.
The path to this file can be set in config_emep.nml, ``femisFile = /my/path/femis.dat``

.. code-block:: text
    :name: base-femis
    :caption: ``femis.dat`` for a base run.

    Name sector  sox  nox  co   voc  nh3  pm25   pmco
    27     0     1.0  1.0  1.0  1.0  1.0  1.0    1.0

This base run example means that there are (1.0), no emission reductions
of |SOx|\ , |NOx|\ , CO, VOC, |NH3|, |PM25| and |PMco| from all sectors in the UK.

-  The first column of the second line represents the country code. (27
   is the code for UK.) The codes for all countries can be found in
   Fortran module ``Country_mod.f90`` or at (http://www.emep.int/grid/country_numbers.txt). 
   The country code must be the same as in the emission files for the given country. Some
   countries and areas are divided into sub-areas in the emission files.
   In this case, one line for each sub-area has to be included into the
   ``femis.dat`` file. Countries and areas where emissions are given for
   sub-areas include the Russian Federation, Germany and all sea areas.
   "0" means all countries.

-  The second column of the second line represents the sector and "0"
   means all sectors. Here one can write the appropriate sector code if
   the emission is reduced only from a specific sector. The description
   of each sector can also be found in the Fortran module ``EmisDef_mod.f90``.

-  The columns under the pollutant names show the emission factors for
   the given pollutants. For example, 0.7 would mean 70% of the original
   emission, thus 30% reduction.

-  The first line is important, it must consist of two words followed
   by the names of the pollutants emitted (in any order).

An example of ``femis.dat`` file describing 50% reduced emission of |NH3|
from sector 10 (the emission from agriculture) in the UK is shown in
:numref:`reduction-femis` and :numref:`reduction-femis_country`.

.. code-block:: text
    :name: reduction-femis
    :caption: ``femis.dat`` for 50% |NH3| reduction from sector 10 over UK.

    Name sector  sox  nox  co   voc  nh3  pm25   pmco
    27    10     1.0  1.0  1.0  1.0  0.5  1.0    1.0

Or alternative syntax using country ISO short names:

.. code-block:: text
    :name: reduction-femis_country
    :caption: ``femis.dat`` for 50% |NH3| reduction from sector 10 over UK.

    Name       sector sox  nox  co   voc  nh3  pm25   pmco
    Country UK  10    1.0  1.0  1.0  1.0  0.5  1.0    1.0

Instead of entire countries, reductions can also be specified by coordinates (and combined with country reductions).
The line with coordinate corrections must start with the keyword ``lonlat``.
The coordinates are given in longitude latitude (min and max and the coordinates of the centre of the gridcells are tested.
Gridcells are either entirely included or entirely reduced, never cut into smaller parts).


.. code-block:: Fortran
    :name: femis
    :caption: ``femis.dat`` example.

    Name                          7  sox  nox  co   voc  nh3  pm25  pmco
    17                            0  1.0  1.0  1.0  1.0  1.0  0.5   0.5
    lonlat 3.3 7.2 50.7 53.5   17 0  1.0  1.0  1.0  1.0  0.0  1.0   1.0


In :numref:`femis`, country with code 17 (NL) will reduce |PM25| and |PM10| emissions by half for all sectors.
Emissions of |NH3| from country with code 17 only, will be removed from the rectangle with longitudes between
3.3 and 7.2 degrees East, and between 50.7 and 53.5 degrees North. Use zero (0) as country code to specify that emissions from all countries should be reduced.



Separate hourly outputs
-----------------------

The ``Base_hour.nc`` and ``Base_uEMEP_hour.nc`` files can become very large. It is possible to split them into one file per day by adding 
the keyword ```HOURLYFILE_ending    = 'JJJ.nc' ```, in the configuration file. ```JJJ```,
will be automatically replaced by the corresponding day of the year (i.e. a number from 1 to 366 giving for instance ``Base_hour_001.nc``).
The full date as in ``Base_hour_20180101.nc`` can be obtained by defining ```HOURLYFILE_ending    = 'YYYYMMDD.nc' ,```. 


.. _`sec-emission-cv-format`:

Using and combining gridded emissions: Country Variable (CV) format
-------------------------------------------------------------------

The gridded emission files are controlled via the ``config_emep.nml`` file. 

In the CV format, emissions are organised in a number of files (Emis_sourceFiles(i_file)), each files containing a number of sources (Emis_sourceFiles(i_file)%source(j_source)).
A source can be any 2D field. It can also be a 3D field, if the third dimension is the sector. The sources are associated to country through a country code or a country code. 
The file must have a ‘lon’ and a ‘lat’ variable, showing longitude and latitudes of each grid point. ‘lon’ and ‘lat’ must be 1D variables if the projection is ‘lon lat’, 2D otherwise.

The code directory contain a Python script that can create an emission file in this format: ``emissions_TXT2ncCV.py``

The file and sources can be characterized by a set of variables. In general these variable can be set by, and in order of increasing priority:
  1. Default value
  2. Global attribute read in the netcdf file
  3. Variable attribute read in the netcdf file
  4. Value set for Emis_sourceFiles(i)%XXX in config_emep.nml
  5. Value set for Emis_sourceFiles(i)%source(s)%XXX in config_emep.nml
  
Exceptions to the priority rule are:
  * maskID cannot be set by attributes in the netcdf file
  * the file and source 'factor' are on top of each other, not replaced
  * boolean parameters (like apply_femis), are used as "and" (i.e. if any is false, the result is false) 

List of file attributes (default in parenthesis):
  - filename (‘NOTSET’) Name of the file (with path)
  - projection (‘lon lat’) Only three categories ‘lon lat’, 'native' or any other (for example ‘Lambert ‘or ‘Stereographic’ would give the same result).
    'native' means that emissions are given in the same grid as the model grid and the data is not interpolated (use it if you can).
  - grid_resolution (an approximate value is computed from the lon and lat, if no value is given).
    It does not need to be exact (cannot be exact on a sphere anyway!).
    This grid_resolution steers the interpolation algorithm;
    A large value will force the code to subdivide each emission gridcell in large number of pieces,
    that are assigned to the model grid. Larger values means smoother interpolation, but more cpu time. 
  - periodicity (‘time’) How often the values are updated.
    Can be ‘yearly’, ‘monthly’, ‘hourly’ or ‘time’. ‘hourly’ or ‘time’ means that the time as defined in the netcdf is used to define when to fetch a new record.
    The timestamp must correspond to the end of the time period of validity.
    For ‘yearly’ monthly timefactors are applied, if a sector is defined.
    For ‘monthly’ and ‘yearly’, an hourly timefactor is applied if a sector is defined. For ‘hourly’ or ‘time’, no additional timefactors are applied. 
  - factor (1.0) multiplicative factor for all sources in the file
  - units ('NOTSET') will be used as default for sources units if set.
  - apply_femis (true) whether to apply the femis reductions to the sources of this file.
  - include_in_local_fractions (true) whether to take sources from this file into account for the local fraction calculations
  - countrycode (-1) will be used as default for sources country code if set. Use rather country_ISO if you can.
  - country_ISO ('NOTSET') will be used as default for sources country ISO code (for example ‘FR’ for France) if set.
  - sector (-1) will be used as default for sources sector if set.
  - SECTORS_NAME or sectorsName ('GNFR_CAMS') . It must be set to one of the SECTORS%name defined ('GNFR_CAMS' for example).
  - mask_ID ('NOTSET') the name of the mask, if you want to apply one.
  - mask_ID_reverse ('NOTSET') the name of the mask, if you want to apply one in the complementary region.
  - country_ISO_incl ('NOTSET') If used, only use emissions from the countries in the list (see example below)
  - country_ISO_excl ('NOTSET') If used, do not take into account emissions from countries in the list (see example below)


List of source attributes:
  - varname (‘NOTSET’) The name as used in the netcdf file
  - species (‘NOTSET’) Either one of the emission group species, as defined in CM_EmisFile.inc (generally sox, nox, pm25, pmco, nh3, co, voc) 
  - factor (1.0) multiplicative factor. Can be used to change units to model definitions. Comes on top of the file multiplicative factors and possibly other factors.
  - units (‘mg/m2/h’) Units *after* the factor multiplication. 
  - countrycode (-1) will be used as default for sources country code if set. Use rather country_ISO if you can.
  - country_ISO (‘N/A’) the country code, as defined in Country_mod.f90 (for example ‘FR’ for France). ‘N/A’ is a valid code, but it does not correspond to any country.
  - apply_femis (true) whether to apply the femis reductions to this source.
  - include_in_local_fractions (true) whether to take this source into account for the local fraction calculations
  - mask_ID ('NOTSET') the name of the mask, if you want to apply one.
  - mask_ID_reverse ('NOTSET') the name of the mask, if you want to apply one in the complementary region.

If a dimension with name "sector" is found, it is assumed that the third dimension of emission variables represent the sector index.

The idea is that only variables that clearly are required in a specific context need to be set;
if the value can be inferred from other information, the code should do it.
Depending of the type of source, not all variables are used.

Note about species: These can be interpreted in one of three categories
  1. emitted species (nox,sox,pm25 ...) with sector (1...11 (or 13)) (“sector species”)
  2. individual species (SO2, NO, NO2, ...) with sector. The species MUST be one of the splitted species.
     These will be treated as one of the “sector species”  from 1. (but not splitted of course).
     Careful with units, it follows the same rules as “sector species”; molecular weight for SO4 for example is considered “as SO2” and NO "as NO2". You can correct the units using the factor attribute.
  3. individual species (SO2, APINENE, O3 ...) without sector (<=0, or not specified).
     No timefactors, vertical realease heights or splits are applied.
     In this case the emissions are summed up in setup_rcemis (not in EmisSet)
     
     Example:

.. code-block:: text
  :caption: Read a field named "emis_tra" from a file "/path/ECLIPSE_V6a_CLE_base_PM25.nc".

  Emis_sourceFiles(1)%filename = '/path/ECLIPSE_V6a_CLE_base_PM25.nc',
  Emis_sourceFiles(1)%sectorsName = 'GNFR',
  Emis_sourceFiles(1)%projection = 'lon lat',
  Emis_sourceFiles(1)%periodicity = 'yearly',  
  Emis_sourceFiles(1)%source(1)%varname='emis_tra',
  Emis_sourceFiles(1)%source(1)%species='pm25',
  Emis_sourceFiles(1)%source(1)%sector=6,
  Emis_sourceFiles(1)%source(1)%factor=1000000.0,!kt->kg
  Emis_sourceFiles(1)%source(1)%units='kg',

If two CV emission files are wished to be combined in a single simulation, for example so that different
emission scenarios can be applied inside and outside of the EMEP region (UNECE excl. NA), one can use
the Emis_sourceFiles(1)%country_ISO_excl(1:) and Emis_sourceFiles(1)%country_ISO_incl(1:) options. 
The below dummy example applies a 2015 baseline scenario inside the EMEP region and 2050 MFR scenario in the rest of the world (ROW),
based on 0.5 degree emission input data files.

.. code-block:: text
  :caption: Include or exclude countries from emission input files.

  Emis_sourceFiles(1)%filename = '/path_to_global_0.5deg_2015_baseline.nc',
  Emis_sourceFiles(2)%filename = '/path_to_global_0.5deg_2050_MFR.nc',
  Emis_sourceFiles(1)%country_ISO_incl(1:48) = "AL", "AM", "AT", "AZ", "BY", "BE", "BA", "BG", "CY", "CZ", "DE",
                                               "DK", "EE", "FI", "FR", "GE", "GR", "HR", "HU", "IS", "IE", "IT",
                                               "LV", "LT", "LU", "MT", "NL", "NO", "PL", "PT", "MD", "RO", "SI", 
                                               "ES", "SE", "CH", "MK", "TR", "SK", "UA", "RS", "ME", "GB", "KG",
                                               "RUSS_EURO", "RUSS_ASIA", "KZT", "FSUA",

  Emis_sourceFiles(2)%country_ISO_excl(1:48) = "AL", "AM", "AT", "AZ", "BY", "BE", "BA", "BG", "CY", "CZ", "DE",
                                               "DK", "EE", "FI", "FR", "GE", "GR", "HR", "HU", "IS", "IE", "IT",
                                               "LV", "LT", "LU", "MT", "NL", "NO", "PL", "PT", "MD", "RO", "SI", 
                                               "ES", "SE", "CH", "MK", "TR", "SK", "UA", "RS", "ME", "GB", "KG",
                                               "RUSS_EURO", "RUSS_ASIA", "KZT", "FSUA",
.. code-block:: Fortran
    :caption: TEST.
    
  Why does this not show?
  Because it is not indented enough!


   
Masks
-----

Typically you have got fine scale emissions for a small region of interest, a city for instance. You may want to remove that area from the coarse scale emissions, and replace it with your own. The mask allows you to define a specific region (the mask).
To define which gridcells to include in your local region, you must find a suitable variable that shows the region of interest. It could be for example the PM emissions in your local area. 

A "mask" can be defined for instance with:

.. code-block:: Fortran
    :caption: Define a cell-fraction mask, example. The mask value represents the fraction of the grid-cell area to be reduced.

    EmisMask(1)%filename = '/mypath/myfile.nc', ! name of the netcdf file to read from
    EmisMask(1)%cdfname  = 'London_PM',  ! name of the variable to read from the file
    EmisMask(1)%type     = 'CELL-FRACTION', ! the mask value is part of the reduction, default
    EmisMask(1)%ID       = 'LONDON',  ! the name you give to that mask
    EmisMask(1)%fac = 0.85, ! muliplicative factor to use. Default 0.0

    ! define species to be read from the emission file separately, as if reading from different files
    Emis_sourceFiles(1:7)%filename=7*'path/Emis_CAMS_v5_1_with_Ref2_0_1_GNFR_CAMS_01005deg_2018_CAMS2_40_U5.nc',
    Emis_sourceFiles(1:7)%species='sox','nox','voc','nh3','co','pm25','pmco',
    ! apply mask to pm25 and pmco emissions
    Emis_sourceFiles(6:7)%mask_ID =2*'LONDON',


.. code-block:: Fortran
    :caption: Define a mask with lower and upper threshold, example

    EmisMask(1)%filename = '/mypath/myfile.nc', ! name of the netcdf file to read from
    EmisMask(1)%cdfname  = 'London_PM',  ! name of the variable to read from the file
    EmisMask(1)%type     = 'THRESHOLD',
    EmisMask(1)%ID       = 'LONDON',  ! the name you give to that mask
    EmisMask(1)%threshold = 1.0E-10, ! the mask is set at any point larger than the threshold
    EmisMask(1)%threshold_max = 100, ! ... and smaller than the threshold_max value (default 1E60)
    EmisMask(1)%fac = 0.85, ! muliplicative factor to use. Default 0.0

    
Several masks can be defined. Each mask is identified by their "ID". If you want to include in the region also the gridcell which are zero, you can set the threshold slightly negative (-1.0E-10), to include the entire region covered by the variable (otherwise zero values would be defined equivalently to outside of region). To include for example all values=23, but not 24, set EmisMask(1)%threshold = 22.5 and EmisMask(1)%threshold_max = 23.5

A mask defines only a region. It is not directly related to any pollutant. 

The masks defined here, will also be applied on files from emis_inputlist (old format), if use_mask is set (but the multiplicative factor is always 0.0 for old format emissions). It is however not possible to set masks by both systems simultaneously. In the old format only one mask can be used at a time. It will be the reunion of all masks produced by the system above (the ID is meaningless and cannot be specified in old format).

To be used with the Local Fractions (see below), one can also define a set of regions defined by integer numbers. For this one must define the ID with the keyword NUMBER: 

.. code-block:: Fortran
    :caption: Define a set of masks with integers, example
    
    EmisMask(1)%filename = '/mypath/myfile.nc', !name of the netcdf file to read from
    EmisMask(1)%cdfname  = 'region_id',  !name of the variable to read from the file. The variable must be an integer!
    EmisMask(1)%type     = 'NUMBER',
    EmisMask(1)%ID       = 'specific-mask-name', 

.. _`sec-nesting`:

Nesting
-------

The model can be run in a large domain and all the concentrations of pollutants stored at fixed intervalls (3 hours typically).
Then we can define a smaller region within the large domain, and rerun the model in the smaller region, using the stored concentrations at the domain boundaries.
It is then possible to make a simulation in a restricted region with fine resolution, but still taking account the effect of pollutants from outside the small region.
This is called nesting. 
The large domain defines the Boundary Conditions (BC, which are only used at the boundaries of the small domain),
and possibly the Initial Conditions (IC, which must be defined everywhere in the small domain, but only for the start date).

Depending on what is wanted different "Nesting modes" can be defined.
The different options are controlled by the ``NEST_MODE_READ`` and ``NEST_MODE_SAVE`` variables in ``Model_config`` in ``config_emep.nml`` file.
The mode options are:

``NEST_MODE_READ``
    'NONE'
        do nothing (default).
    'START'
        read at the start of run.
    'RESTART'
        read at the start of run, but do not overwrite BC. Use this option when you restart a run
without changing the domain size ("checkpoint-restart"). The difference with "START" is only that BC are not overwritten. 
(This is important in the case of restarting exactly at the start of a month (at time 00:00, the 1st),
because the nest file contains the BC from preceding month, and if 'START' is used it will use those for the new month, 
i.e. overwrite new BC with old ones.)
    ‘NHOUR’
        read at given ``NEST_NHOURREAD`` hourly intervals, if the file is found.
        ``NEST_NHOURREAD`` is set in ``Model_config`` and should be an integer fraction of 24.
    ‘MONTH’
        read at start of each month.
``NEST_MODE_SAVE``
    'NONE'
        do nothing (default).
    'END'
        write at end of run.
    'NHOUR'
        write at given `NEST_NHOURSAVE` hourly intervals.
        ``NEST_NHOURSAVE`` is set in ``Model_config`` and should be an integer fraction of 24.
    ‘MONTH’
        write after each month (used for checkpoint/restart for instance).

The name of the file to write to is defined by ``NEST_template_write`` (also in ``Model_config``).
The name of the file to read to IC data from is defined by ``NEST_template_read_3D``.
The name of the file to read to BC data from is defined by ``NEST_template_read_BC``; it can be the same file as ``NEST_template_read_3D``.
If for example NEST_NHOURSAVE=3, but  NEST_NHOURREAD = 1, the data is interpolated in time to get a smooth transition between the 3-hourly values (recommended).

In addition, 3D outputs can be asked for at specific date. Example:

.. code-block:: text
    :name: nest-dump-config
    :caption: Write 3D at specific dates example.

    &Model_config
    [...]
      NEST_OUTDATE_NDUMP    = 2,          ! choose two dates
      NEST_outdate(1)= 2019,01,01,00,     ! save at end of the 2018 year
      NEST_outdate(2)= 2018,05,17,00,     ! save also 17th of May 2018
    &end
 

Example write BCs
~~~~~~~~~~~~~~~~~

:numref:`nest-write-config` shows an example to write every 3 hours into
daily Nest/BC files. Output file name is defined by ``NEST_template_write`` ('BC_YYYYMMDD.nc'),
where ‘YYYYMMDD’ is replaced by corresponding date.

All advected model variables will be written out for the sub-domain defined
by ``NEST_out_DOMAIN`` (\ :math:`x=60,\ldots,107; y=11,\ldots,58`\ ).
If no ``NEST_out_DOMAIN`` is given, the entire model rundomain will be written out.

.. code-block:: text
    :name: nest-write-config
    :caption: Write BCs configuration example.

    &Model_config
    [...]
      NEST_MODE_READ        = 'NONE',          ! do not read BC
      NEST_MODE_SAVE        = 'NHOUR',         ! write BCs
      NEST_NHOURSAVE        = 3,               !  every 3 hours
      NEST_template_write   = 'BC_YYYYMMDD.nc' !  to your (daily) BC output file
                                          !  (8 records per file)
      !-------- Sub domain for write modes
      NEST_out_DOMAIN       = 60,107,11,58,    ! istart,iend,jstart,jend
    &end
    

Read BCs produced by a previous EMEP MSC-W model run
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:numref:`nest-read-config` shows an example to read every 1 hours from
the Nest/BC files created previously by running :numref:`nest-write-config`.
Please note that the model sub-domain for a nested run is set by ``RUNDOMAIN``,
as shown in :numref:`config-emep`.

.. code-block:: text
    :name: nest-read-config
    :caption: Read BCs configuration example.

    &Model_config
    [...]
      NEST_MODE_READ        = 'NHOUR',         ! read external BC
      NEST_NHOURREAD        = 1,               !   every hour
      NEST_template_read_BC = 'BC_YYYYMMDD.nc' !   your (daily) BC input file
      NEST_MODE_SAVE        = 'NONE',          ! do not write BCs
    &end


Note that ``NEST_NHOURREAD`` can (and should) be smaller than the value used for ``NEST_NHOURSAVE``.
The values between saved dates will then be interpolated in time, giving a smoother transition.


Reduce the size of BC files
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The size of the files obtained in a nesting configuration can be very large if the out_DOMAIN is large.
If the inner domain is known in advance, only the part matching exactly the part needed to construct the BC of the small domain can be stored.
Define ``NEST_MET_inner`` in ``&Model_config``, which should be a link to any metdata of the inner grid;
it will only be used to define the projection parameters of the inner grid (i.e. dates and other content do not matter).

.. code-block:: text
    :name: nest-write-inner
    :caption: Inner domain options for nested BC output example.

    &Model_config
    [...]
      NEST_MET_inner = 'inner_domain/wrfout_d03_2015-01-01_00:00:00',
      NEST_RUNDOMAIN_inner = 12, 136, 100, 300,
    &end
        
You cannot use the implicit dates ("YYYY" etc.); you must put explicit numbers. 
Note that the file will have the same dimensions, but zeros are put into the unused parts.
The NetCDF internal compression will take care of reducing the actual size, as measured by used disc space.

If a BC file has been created using the ``NEST_MET_inner`` method, it cannot be used for initializating concentrations at the start of the run.
A separate file has to been created for initializations.
This file can then be used by the inner grid by defining ``NEST_template_read_3D`` in ``config_emep.nml``.


Read external BCs
~~~~~~~~~~~~~~~~~

So far only BC created by the model itself have been used.Reading external BCs,
i.e. produced by other means (another model for example) is more involved.
The chemical species may be different and the vertical levels also.
The vertical axis and variables in the file need then to be mapped to the corresponding model variables. 

:numref:`nest-mybc-config` shows an example to read every 3 hours from an external
BC file. The model will read 3 variables from ``MyBC.nc``: |O3|, NO, and |NO2|.
The maping between the ``MyBC.nc`` variables and the corresponding model variables
is defined in the ``ExternalBICs_bc`` namelist.

.. code-block:: text
    :name: nest-mybc-config
    :caption: External BCs configuration example.

    &Model_config
    [...]
      NEST_MODE_READ        = 'NHOUR',    ! read external BC
      NEST_NHOURREAD        = 3,          !   every 3 hours
      NEST_template_read_BC = 'MyBC.nc'   !   from your BC input file
      NEST_MODE_SAVE        = 'NONE',     ! do not write BCs
      USE_EXTERNAL_BIC = T,
      EXTERNAL_BIC_NAME= 'MyBICScenario', ! variable mapping, see ExternalBICs_bc
      TOP_BC           = T,               ! use top BC also from your BC file
      filename_eta     = 'filename.zaxis',! text file containing the
                                          !   description of your BC file
    &end
    &ExternalBICs_bc
      description='MyBICScenario','Version name',3,  ! name,version,size
      map_bc=! emep,external,frac,wanted,found,IXADV,
        'O3' ,'O3_VMR_inst' ,1.0,T,F,-1,
        'NO' ,'NO_VMR_inst' ,1.0,T,F,-1,
        'NO2','NO2_VMR_inst',1.0,T,F,-1,
    &end


Vertical coordinate
___________________

In order the determine the vertical levels on the external BC file
('MyBC.nc' in :numref:`nest-mybc-config`),
the following checks will take place in the following order:

#. :math:`\eta` (eta) coordinate:
       If the variable :math:`hyam` (hybrid ``a`` coefficient at layer midpoint) is found,
       :math:`\eta` is calculated from :math:`hyam` and :math:`hyam` variables on the file.
#. :math:`\sigma` (sigma) coordinate:
       If the vertical level is indexed by variable ``k``.
#. :math:`\eta` coordinate:
       If the file defined by variable ``filename_eta`` exist
       ('filename.zaxis' in :numref:`nest-mybc-config`),
       :math:`\eta` is derived from ``vct`` and ``vctsize`` variables on the file.
#. pressure coordinate:
       If the vertical level is indexed by variable ``lev``.

Independent of the coordinates of the BC file,
the BC levels will be interpolated into EMEP model levels.
If the BC file level structure is not recognized,
and there is no ``filename_eta`` provided,
the model will write an error message and stop execution.

An example of the ``filename_eta`` for EMEP model levels is given below.
Here the ``vct`` variable describes the model level boundaries in hybrid eta coordinate:

.. literalinclude:: filename.zaxis
    :caption: ``filename_eta`` for EMEP model levels.

``vct`` is the vertical coordinate table describing the hybrid ``a`` and ``b``
paramters at the layer interfaces in :math:`\eta` coordinate system
(:math:`hyai` and :math:`hybi`). They must respect the following constraints:

* :math:`hyai_1=0` and  :math:`hybi_1=1`;
* :math:`hyai_0=P_t` and :math:`hybi_0=0`.
* :math:`P_t` is the pressure at top of the model domain.

In this file, the first 21 values in ``vct`` represent :math:`hyai`
and the remaining 21 represent :math:`hybi` values in hPa.

Variable mapping
________________

The variables to be used from the external boundary condition data are
given in the ``ExternalBICS_bc`` namelist in the ``config_emep.nml`` file.
In :numref:`nest-mybc-config`, `map_bc` maps 3 variables in the 'MyBC.nc' file to
3 model variables (|O3|, NO, and |NO2|\ ).
On each line of `map_bc` contains the 6 elements:

#. Variable name in EMEP MSC-W model, e.g. `'O3'`.
#. Variable name in the External BC data file, e.g. `'O3_VMR_inst'`.
#. External BC component to EMEP component fraction, e.g. `1.0`.
#. Is this component wanted or not, set to `T` or `.true.` to read the variable.
#. Was the BC variable found on the file, will be set by the model on run time.
#. Index of the advected model variable, will be set by the model on run time.

The fraction is helpful, when one has to map a variable that is explicitly not in the model
but a fraction of that variable can be mapped to a matching variable in the model.

Unit conversions are delegated to the ``Units_mod.f90`` module.
The supported units are: |ugm3|\ , |ugSm3|\ , |ugNm3|\ , |ugCm3|\ , ppb,
mixing ratio (mol/mol) and mass mixing ratio (kg/kg).

If the BC data has different units,
either convert them into one of the above mentioned units in pre-processing
or add the respective conversion factor in the module ``Units_mod.f90``.


config: Europe or Generic/Global?
-------------------------

The EMEP model has traditionally been run on the EMEP grid covering
Europe, and using meteorology from the ECMWF IFS model. In this environment,
we typically set several configuration variables to make use of
Euro-specific data. In other regions it is better to make use of the
model's generic ('global') settings, which will ensure better handling of
vegetation (eg LAI) changes, convection, and various emission
settings. For the EMEP grid, the relevant USES config flag settings are

.. code-block:: text
  :caption: Typical European/EMEP settings.

  USES%DEGREEDAY_FACTORS = T,    ! though F is okay too
  USES%PFT_MAPS = F,
  USES%CONVECTION = F, 
  USES%ROADDUST = F,
  Vertical_levelsFile = 'DataDir/Vertical_levels20_EC.txt'

while for the generic domain settings these are

.. code-block:: text
  :caption: Typical non-European/global settings.

  USES%DEGREEDAY_FACTORS = F
  USES%PFT_MAPS = T,    ! PFT LAI tests
  USES%CONVECTION = T,
  USES%ROADDUST = F,
  Vertical_levelsFile = 'DataDir/Vertical_levels19_EC.txt',
  timefacs%Monthly = 'GRIDDED',

While the meaning of the above flags is described in more detail below, since EMEP model version 
rv5.4 the sets of above flags are controlled by the USES%DOMAIN_SETUP_TYPE flag, and no longer included
as separate USES flags in the config namelist. The DOMAIN_SETUP_TYPE flag currently
has three options: 'EMEPDOMAIN', 'GenericDOMAIN', and 'Custom_Config'. The 'EMEPDOMAIN' and 'GenericDOMAIN' flags
trigger their respective USES options described above, and overwrite any USES settings for these flags
included in the config namelist (by default, these flags are therefore not included in the namelist). If a user wishes 
to use custom domain settings, using the 'Custom_Config' domain setting re-activates the use of the above USES
flags in the config namelist (i.e., these are not overwritten with what would be specified through the 'EMEPDOMAIN'
and 'GenericDOMAIN' flags). Alternatively, one could for example use 'EMEPDOMAIN', while also specifying
EMEP_DOMAIN_SETUP%USES_PFTMAPS = T in the config namelist, if the plan is to use PFT maps for European scale modelling. 
Additionally, one may specify EMEP_DOMAIN_SETUP%USES_DEGREEDAYS = F to run the simulation without degreedays.
To avoid confusion as to which domain settings are applied, users are forced to choose either one of the three
above domain settings. 

The ``DEGREEDAY_FACTORS`` setting triggers the use of degree-days in
controlling residential combustion (GNFR C/SNAP2) emissions. This requires pre-processed files of heating degree days.
Such files (DegreeDayFactors.nc) can be produced from any meteorology, but the difference in
results even for Europe is not too significant. In other regions of the
world emissions from residential combustion may not be as dependent on degree-days as in
Europe, and so this setting should probably be false.

``PFT_MAPS=T`` triggers the use of a global file which provides monthly
variations in leaf area index (LAI) for different vegetation types. This
controls deposition and biogenic VOC emission parameters. In Europe, a
simpler latitude-dependent system is used (based upon DO3SE), and so
``PFT_MAPS`` should be set F.

``CONVECTION`` is difficult. In principle, all models runs should use the T
setting, but for Europe we find it degrades the model results too much
and we use F. The problem is likely that the sub-grid processes behind
convection are so complex and the paramererisation is very uncertain.
Note also that in Config_module we have the default setting ``CONVECTION_FACTOR=0.33``,
which may be changed to allow more or less influence of this variable.

``ROADDUST`` includes dust produced by road traffic based on input files broadly
covering the EMEP modelling domain (see Input section). However, this input is now deprecated
and both the Generic and EMEP domain setups set this flag to False.

``Vertical_levelsFile`` controls the number of verticl levels, with the Generic domain 
setting employing 19 vertical levels rather than 20 for the EMEP domain. For regions outside 
of Europe high orography and tall trees can make it more scientifically sound to use a 
thicker surface layer, with the bottom two layers of the 20-level file being combined
into a single layer being approximately 90 m thick.

``timefacs%Monthly`` Outside of Europe, it is also typical to use time-factors derived
from the global CAMS-TEMPO data set.

config: SOILNOX
---------------


.. code-block:: text
  :caption: Soil-NO settings.

  USES%SOILNOX_METHOD = 'NoFert',  ! If using ECLIPSEv6 or EMEP European
  USES%SOILNOX_METHOD = 'Total',   ! If using ECLIPSEv5

By default the model makes use of global 0.5 degree data from the CAMS2-61 project (cf Simpson et al., Ch.8 in
Denier van der Gon, 2023, doi:10.24380/q2si-ti6i,
https://atmosphere.copernicus.eu/node/1054), but the user needs
to specify the data to be used from this system. 
The choice, between ``Total``  and ``NoFert`` depends
on the anthropogenic emission inventory in use.

*IMPORTANT*:

If this anthropogenic inventory
already includes fertlizer-induced soil NO emissions, then the  ``NoFert`` data from
the soil-NO system is used. This case applies to EMEP/CEIP emissions,  HTAPv3, EDGAR and
CAMS-GLOB-ANT.
The latest ECLIPSE inventories (v6+) also incude fertlizer-induced emssions

If fertlizer-induced soil NO emissions are not included in the inventory,
then choose ``Total``. This case applies to CAMS-REG-ANT and some earlier ECLIPSE (v5
and earlier) inventories.


If the earlier soil-NOx methods from Simpson et al. (2012) are prefered, then these
can be set using:

.. code-block:: text
  :caption: Alternative (ACP 2012) Soil-NO settings, for European runs.

  USES%SOILNOX_METHOD = 'ACP2012',  !  Europe only
  NdepFile            = 'DataDir/AnnualNdep_PS50x_EECCA2005_2009.nc',









Other less used options
-----------------------

There are many internal settings that are set to a default value by the model.
These default values can often be overridden by setting specific values for keywords in ``config_emep.nml``.

The ``RUNDOMAIN`` is divided by the model into subdomains which are assigned to each processor.
Normally this partioning is done such that the X and Y direction are divided into approximately equal number of parts.
For runs in lat lon projection containing poles the Y division is done into one or two parts,
so that each processor has the same share of pole regions.
The default partitioning can be overrided using the ``DOMAIN_DECOM_MODE`` parameter in ``config_emep.nml``.
Recognized values are: 
`'X*Y'`, `'XY'`, `'X*1'`, `'Y=1'`, `'X'`, `'1*Y'`, `'X=1'`, `'Y'`, `'2*Y'`, `'X=2'`, `'X*2'`, `'Y=2'`.
See also in ``Par_mod.f90`` for details.

The main timestep parameter ``dt_advec`` can be set manually in config_emep.nml (in seconds, real number).
It must be a factor of 3600.

``JUMPOVER29FEB`` : if set to T , will not treat the 29th of February even for leap years.
(Useful if that date is missing, for instance in Climate runs).

``NETCDF_DEFLATE_LEVEL``: The level of compression used in the NetCDF output files (integer).
Negative values means netcdf3 format.

``END_OF_EMEPDAY``: (integer) Hour when to start and end a "day". For daily averages gives for example averages from 06:00 to 06:00 next day, if ``END_OF_EMEPDAY = 6``.

``USES%SKIP_INCOMPLETE_OUTPUT`` : if T, will skip daily/montly/fullrun output for runs under 1/28/181 days.

.. _`sec-emission-own-sectors`:

Defining own sectors
--------------------

Emissions can be assigned to a sector. A sector defines three properties:


1. Split into species
2. Emission height release distribution
3. Timefactors

The set of splits is defined in files "emissplit.defaults.POLL" and "emissplit.specials.POLL".
The release height distribution can be chosen among the 6 predefined values, or defined in config_emep.nml.
The timefactors are defined in "MonthlyFac.POLL" "DailyFac.POLL" and "HourlyFacc.INERIS".

Which height/split/timefac is chosen for a given sector can be also be controlled by the user.
GNFR_CAMS sectors are predefined (SNAP also, but we recommend to use GNFR_CAMS).
The values for split, emission release height and timefactors can be defined through settings in the config_emep.nml settings. The following will reproduce the default settings:

.. code-block:: Fortran
    :name: gnfr-cams-sectors
    :caption: Sector definitions equivalent to the predefined sectors for GNFR_CAMS.

    SECTORS_ADD(1) = 'GNFR_CAMS', 'GNFR_A',  'sec01',  1, 1,  1, 'Public Power', 'ALL',
    SECTORS_ADD(2) = 'GNFR_CAMS', 'GNFR_B',  'sec02',  3, 3,  2, 'Industry', 'ALL',
    SECTORS_ADD(3) = 'GNFR_CAMS', 'GNFR_C',  'sec03',  2, 2,  3, 'OtherStationaryComb', 'ALL',
    SECTORS_ADD(4) = 'GNFR_CAMS', 'GNFR_D',  'sec04',  5, 5,  4, 'Fugitive', 'ALL',
    SECTORS_ADD(5) = 'GNFR_CAMS', 'GNFR_E',  'sec05',  6, 2,  5, 'Solvents', 'ALL',
    SECTORS_ADD(6) = 'GNFR_CAMS', 'GNFR_F',  'sec06',  7, 2,  6, 'RoadTransport', 'ALL',
    SECTORS_ADD(7) = 'GNFR_CAMS', 'GNFR_G',  'sec07',  8, 2,  7, 'Shipping', 'ALL',
    SECTORS_ADD(8) = 'GNFR_CAMS', 'GNFR_H',  'sec08',  8, 7,  8, 'Aviation', 'ALL',
    SECTORS_ADD(9) = 'GNFR_CAMS', 'GNFR_I',  'sec09',  8, 2,  9, 'Offroad', 'ALL',
    SECTORS_ADD(10) = 'GNFR_CAMS', 'GNFR_J', 'sec10',  9, 6, 10, 'Waste', 'ALL',
    SECTORS_ADD(11) = 'GNFR_CAMS', 'GNFR_K', 'sec11', 10, 2, 11, 'AgriLivestock', 'ALL',
    SECTORS_ADD(12) = 'GNFR_CAMS', 'GNFR_L', 'sec12', 10, 2, 12, 'AgriOther', 'ALL',
    SECTORS_ADD(13) = 'GNFR_CAMS', 'GNFR_M', 'sec13',  5, 5, 13, 'Other', 'ALL',
    SECTORS_ADD(14) = 'GNFR_CAMS', 'GNFR_A1','sec14',  1, 1,  1, 'PublicPower_Point', 'ALL',
    SECTORS_ADD(15) = 'GNFR_CAMS', 'GNFR_A2','sec15',  1, 3,  1, 'PublicPower_Area', 'ALL',
    SECTORS_ADD(16) = 'GNFR_CAMS', 'GNFR_F1','sec16',  7, 2, 16, 'RoadTransportExhaustGasoline', 'ALL',
    SECTORS_ADD(17) = 'GNFR_CAMS', 'GNFR_F2','sec17',  7, 2, 17, 'RoadTransportExhaustDiesel', 'ALL',
    SECTORS_ADD(18) = 'GNFR_CAMS', 'GNFR_F3','sec18',  7, 2, 18, 'RoadTransportExhaustLPGgas', 'ALL',
    SECTORS_ADD(19) = 'GNFR_CAMS', 'GNFR_F4','sec19',  7, 2, 19, 'RoadTransportNonExhaustOther', 'ALL',


The name in the first column (GNFR_CAMS) should match the sector name defined in the emission file, or be given in the config_emep.nml (see example below).
The second name can be chosen by the user (it will be used if SecEmisOutWanted is T)), the third name is the variable name ending in the netcdf emission file (the start being the name of the pollutant) only implemented for "fraction" format, the first number refers to the index used for time factors, the second number the index in the Emis_h array, the third number is the index used in the split files. The long name in the second last column is a longer description, and the last column refers to the species to be included for this sector (if they exist in the emission file).

From version 4.45 and onwards, sector type or names are compulsory for each emission file. This can be given either in the config_emep.nml or in the NetCDF file.

Here is an example of how to define a new sector with a new height distribution, used by emissions given in a separate file.


.. code-block:: Fortran
    :caption: Settings for defining the pm 2.5 emissions from the file MyEmis.nc with variable name 'pm25_MyCar', with emission released between 20 and 50 meters.

    Emis_Zlevels(1:)20.0,   50.0,   92.0,  184.0,  324.0,  522.0,  781.0, 1106.0,
    Emis_h(1:,1) = 0.000,  0.000,  0.000,  0.003,  0.147,  0.400,  0.300,  0.150, 
    Emis_h(1:,2) = 1.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
    Emis_h(1:,3) = 0.060,  0.067,  0.093,  0.750,  0.030,  0.000,  0.000,  0.000, 
    Emis_h(1:,4) = 0.050,  0.063,  0.087,  0.700,  0.100,  0.000,  0.000,  0.000, 
    Emis_h(1:,5) = 0.020,  0.034,  0.046,  0.600,  0.300,  0.000,  0.000,  0.000,
    Emis_h(1:,6) = 0.000,  0.000,  0.000,  0.410,  0.570,  0.020,  0.000,  0.000,
    Emis_h(1:,7) = 0.200,  0.300,  0.020,  0.044,  0.066,  0.094,  0.123,  0.153, 
    Emis_h(1:,8) = 0.000,  1.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, 
      
    emis_inputlist(1)%name='GNFR.nc',
    emis_inputlist(2)%name='MyEmis.nc',
    emis_inputlist(2)%sector='MyNewSector',
    SECTORS_ADD(1) = 'MyNewSector', 'MyTestSector',  'MyCar',  7, 8,  6, 'Special car exhaust', 'pm25',
 
Note that if you define new splits, you must include defaults values in all the default files (even if they are overwritten by the specials).


Local Fractions
---------------

The Local Fraction method allows to give information about where the pollutants come from. There are two main "branches", tracking of primary particles, and "generalized LF", which also includes chemical transformations between species. The values for the local fractions will be outputted in separate files (file names including  "_LF_"). 




Local Fractions for primary particles
-------------------------------------

Primary particles are simplest; we can imagine that pollutants from each source are tracked separately (it is what we do in the code). What is included as a source is defined by the settings in config_emep.nml.

Mainly two types of sources:

* grid to grid (aka relative), which considers each individual gridcell separately. To keep the amount of data reasonable, the pollutants are tracked only up to a given distance (10 gridcells in each direction for example). 
* countries, possibly limited to a specific sector. (SR type runs)

Below an example of grid to grid LF which are the type used by the urban EMEP, uEMEP.

config_emep.nml settings:

.. code-block:: Fortran
    :caption: Local Fractions flag example

    USES%LocalFractions = T, ! T for computing Local Fractions, F otherwise
    !Local Fractions frequency of output (separate file for each). Can be any of: YEAR, MONTH, DAY, HOUR, HOUR_INST 
    lf_set%YEAR = T, !average value for full run in output
    lf_set%Nvert = 14, !How many vertical level to include in treatment. Should be higher than highest emissions

    !NB: for now lf_src(1)%dist will be used for all sources
    lf_src(1)%dist = 5,  !how far the neighbors can be in each direction (NB: high cost for large dist)
    
    !Local Fractions pollutants and sectors to include:
    lf_src(1)%species="pm25", ! any of EMIS_File: "sox ", "nox ", "co  ", "voc ", "nh3 ", "pm25", "pmco"
    lf_src(1)%sector=0, !0 means sum of all sectors
    lf_src(1)%make_fracsum=T,
    lf_src(2)%species="pm25", ! any of EMIS_File: "sox ", "nox ", "co  ", "voc ", "nh3 ", "pm25", "pmco"
    lf_src(2)%sector=0, !0 means sum of all sectors
    lf_src(2)%res=5, ! sources are squares with size 5x5 gridcells
    lf_src(2)%make_fracsum=T,
    lf_src(3)%species="nox ",
    lf_src(3)%sector=0,
    lf_src(4)%species="nox ",
    lf_src(4)%sector=8,
    
    ! lf_set%DOMAIN = 370, 420, 270, 320, !which domain to include in output. Will save disk, but not CPU to reduce.


If one wants to include many species, sectors and res values, without writing one entry per source, one can use the following syntax:

.. code-block:: Fortran
    :caption: Local Fractions sectors arrays example

    USES%LocalFractions = T, ! T for computing Local Fractions
    lf_species(1:2)%name = 'pm25','nox',
    lf_species(1)%sectors(1:) = 0, 1, 2, 8,
    lf_species(1)%res(1:) = 1, 4,
    lf_species(2)%sectors(1:) = 0, 1, 2, 8,
    lf_species(2)%res(1:) = 1, 4,

The corresponding lf_src values will then be added to the already defined lf_src (2*4*2 = 16 new sources in this example). You can not with this syntax combine different res for different sectors for the same species.

NB: by default, the GNFR sector 6 also includes 16,17,18 and 19; the GNFR sector 1 also includes 14 and 15.

Note that the files can be very large if hourly outputs and/or many neighbors are requested.



Local fractions can also be used to make traditional Source Receptor (or blame) matrices, in a single run. (the %res flag is then not used)

.. code-block:: Fortran
    :caption: Local Fractions Country source receptor type example

    USES%LocalFractions = T, ! T for computing Local Fractions
    !Local Fractions frequency of output (separate file for each). Can be any of: YEAR, MONTH, DAY, HOUR, HOUR_INST 
    lf_set%YEAR = T, !average value for full run in output 
    lf_set%Nvert = 14, !How many vertical level to include in treatment. Should be higher than highest emissions
    
    !Local Fractions pollutants and sectors to include:
    lf_src(1)%species="pm25", ! any of EMIS_File: "sox ", "nox ", "co  ", "voc ", "nh3 ", "pm25", "pmco"
    lf_src(1)%type='country',  ! Means make country style SR
    lf_src(1)%drydep=T, !means make country dry deposition maps too
    lf_src(1)%wetdep=T, !means make country wet deposition maps too

    lf_src(2)%species="pmco", ! any of EMIS_File: "sox ", "nox ", "co  ", "voc ", "nh3 ", "pm25", "pmco"
    lf_src(2)%type='country',  ! Means make country style SR
    lf_src(2)%drydep=T, !means make country dry deposition maps too
    lf_src(2)%wetdep=T, !means make country wet deposition maps too

    ! Specify which countries and sectors
    lf_country%sector_list(1:8)=0,1,2,3,4,5,6,7,
    lf_country%list(1:20)='FR','IT','DE','ES','NO','NL','SE','PL','AT','BE','BG','DK','FI','GR','HU','PT','RO','CH','TR','GB',
    lf_country%group(1)%name='NORDIC', !any name given to the group (used as output name)
    lf_country%group(1)%list(1:)='NO','DK','SE','FI', ! countries included in the group
    
Instead of defining countries in the emission files, one can define "source regions" in a separate netcdf file. Each region must have an integer value. The values to be included as a source region are then specified by a list or as the minimum and maximum value to include (defining both a list and an intervall is allowed, but only one list and one intervall). The maskfile should be defined too. For example to include masks 2,5,8,12 and all mask between 100 and 357 (included):

.. code-block:: Fortran
    :caption: Local Fractions mask regions source receptor example

    EmisMask(1)%filename='municip_mask/municip_mask_500m.nc',
    EmisMask(1)%cdfname =region_id',
    EmisMask(1)%type    ='NUMBER',
    EmisMask(1)%ID      ='municip_mask',

    lf_country%mask_val(1:) = 2,5,8,12, 
    lf_country%mask_val_min = 100,
    lf_country%mask_val_max = 357,
    lf_country%sector_list(1:1)=0, ! sum of all sectors
    lf_country%cellmask_name(1:1) = 'municip_mask',
    lf_src(1)%species="pm25",
    lf_src(1)%type='country',
    mask2name(1)='Oslo',     !map id numbers to names (for use in output)
    mask2name(2)='Eigersund',
    
If a value is within the min and max range, but does not appear in the mask file, it will not be taken into account (meaning it is ok to specify a range that covers all the masks values, even if some values are not defined on the mask) 
    
use mask file with fraction of emissions in each gridcell: 

.. code-block:: Fortran
    :caption: Local Fractions using mask with fractions of gridcells

    EmisMask(1)%filename = '/ec/res4/hpcperm/fan/Data/Masks/cameo_city_masks.nc',
    EmisMask(1)%cdfname = 'TouHan',
    EmisMask(1)%type = 'CELL-FRACTION',
    EmisMask(1)%ID = 'ToHa',

    EmisMask(2)%filename = '/ec/res4/hpcperm/fan/Data/Masks/cameo_city_masks.nc',
    EmisMask(2)%cdfname = 'PorUtr',
    EmisMask(2)%type = 'CELL-FRACTION',
    EmisMask(2)%ID = 'PoUt',

    lf_country%cellmask_name(1) =  'ToHa',

In this example two masks are defined (those can be used for traditional SR runs too), and only the mask with name "PoUt" is used as a "country" in the LF run.


Local Fractions for tracing primary natural emissions (under development)
-------------------------------------------------------------------------
Natural emissions are assumed emitted from surface (in current version). 
So far only the DMS emissions can be tracked. 
To include DMS outputs:

.. code-block:: Fortran
    :caption: Example for DMS output

    lf_src(1)%name = 'DMS',
    lf_src(1)%dist = 2, !will track over up to 2 gridcells in all directions
    lf_src(1)%nhour = 1, ! will track separately emissions every 1 hour, and reset every 24 hours
    lf_set%HOUR_INST = T, ! output instantaneous values every hour

Local Fractions for Sensibilities with full chemistry
-----------------------------------------------------

When species from different sources mixes, it is not trivial to define the source of the resulting species. Instead of "fractions" we define sensibilities (or generalized Local Fractions); sensibilities tells how much a species would vary if a given source would change slightly. Mathematically well defined as a derivative of species wrt to emission source. We define the units such that the value of the sensibility gives the change of concentration extrapolated to a 100% emission change. For linear or primary particles this means that the concentration are simply the amount of pollutant coming from that source (i.e. the local fraction multiplied by the base concentration).


The cpu cost is high, approximatively 20 times the cost without this option (independently of the number of sources tracked, up to 100 or so sources). To use this option the fortran code must be prepared with the script ``utils/mk.LF_Chem``, or overwrite ``CM_Reactions1.inc`` with ``CM_Reactions1_LF.inc`` and ``CM_Reactions2.inc`` with ``CM_Reactions2_LF.inc``. Example of config settings:

.. code-block:: Fortran
    :caption: Local Fractions Country source receptor type example

    USES%LocalFractions = T, ! T for computing Local Fractions
    lf_set%YEAR = T, !average value for full run in output 
    lf_set%Nvert = 14, !How many vertical level to include in treatment. Should be higher than highest emissions
    lf_set%full_chem = T, ! to indicate that all species must be included
    lf_set%Nfullchem_emis = 4, ! can be 1,2 or 4. 1 reduces nox, voc, sox, nh3 together; 2 reduces nox, voc separately; 4 reduces nox, voc, sox, nh3 separately

    ! Specify which countries and sectors
    lf_country%sector_list(1:)=0,1,8,
    lf_sector_groups(1)%name='Low', !any name you want to give to this group
    lf_sector_groups(1)%list(1:)=3,5,6,9,11,12,16,17,18,19, !the sector indices included in this group
    lf_country%list(1:20)='FR','IT','DE','ES','NO','NL','SE','PL','AT','BE','BG','DK','FI','GR','HU','PT','RO','CH','TR','GB',
    lf_country%group(1)%name='NORDIC', !any name given to the group (used as output name)
    lf_country%group(1)%list(1:)='NO','DK','SE','FI', ! countries included in the group


    ! Specify which species or group of species must be outputted
    lf_spec_out(1)%name='O3', !single species
    lf_spec_out(2)%name='NO', !single species
    lf_spec_out(3)%name='NO2', !single species
    lf_spec_out(4)%name='SIA',!group of species
    lf_spec_out(4)%species(1:)='SO4','NO3_f','NO3_c','NH4_f', !species to include in the group
    lf_spec_out(4)%species_fac(1:)=1.0, 1.0, 1.0, 1.0,        !weights (default 1.0)
    lf_spec_out(5)%name='NOx',!group of species
    lf_spec_out(5)%species(1:)='NO','NO2', !species to include in the group
    lf_spec_out(6)%name='PM25_rh50',
    lf_spec_out(6)%species(1:)='NO3_c','PM_WATER','SO4','NO3_f','NH4_f','ASOC_ng1e2','ASOC_ug1', 'ASOC_ug10', 'ASOC_ug1e2', 'ASOC_ug1e3','non_C_ASOA_ng1e2', 'non_C_ASOA_ug1', 'non_C_ASOA_ug10','non_C_ASOA_ug1e2', 'non_C_ASOA_ug1e3',
    lf_spec_out(6)%species_fac(1)=0.202764826447516,
    lf_spec_out(7)%name='DDEP_OXN',
    lf_spec_out(7)%species(1:)='NO2', 'N2O5', 'HONO', 'HNO3', 'HO2NO2', 'SC4H9NO3', 'NALD', 'ISON', 'PAN', 'MPAN', 'NO3_f', 'NO3_c', 'shipNOx',
    lf_spec_out(7)%DryDep=T,
    lf_spec_out(8)%name='WDEP_SOX',
    lf_spec_out(8)%species(1:)='SO2','SO4',
    lf_spec_out(8)%WetDep=T,
    lf_spec_out(9)%name='POD1_IAM_DF',
    lf_spec_out(9)%species(1:)='POD',
    lf_spec_out(9)%DryDep=T,
    lf_spec_out(10)%name='pm25',
    lf_spec_out(11)%name='PM_WATER',

    lf_set%MDA8 = T, !special: make AvgMDA8_6month and SOMO35 NB: requires that O3 is outputted too
    

Technical
=========

Time steps
----------
There are three types of timesteps in the model: ``METSTEP``, ``dt_advec`` and ``dtchem``.

``METSTEP`` is the time interval between two readings of the meteorological file. It is not set by the user, but uses whatever is used as meteorological input. It must be a divisor of 24 hours. Three hours is standard, but other intervals are used too.

``dt_advec`` is the "time splitting" interval, i.e. the time between two sequences of advection-chemsitry/emissions-deposition. (Note that the name is misleading, the advection can have a smaller internal timestep if it is requird by the Courant number). The value of dt_advec must be an entire fraction of one hour (units seconds). It is determined by the model using the grid resolution. For a grid resolution larger than 61 km it is 1800s, or 1200s if larger than 21km, or 900s if larger than 11km, or 600s if larger than 6km, or 300s if larger than 2km or 100s if smaller than 2km. Those value can be overriden by defining it in ``config_emep.nml``. You can check the value in the standard output ("advection time step (dt_advec) set to: 1200 seconds").


The values of ``dtchem`` are the internal Chemistry timesteps. It is variable: usually the first 5 steps last 20 seconds, then followed by 10 larger timesteps, so that the sum of all ``dtchem`` timesteps is exactly equal to ``dt_advec``. They are printed out in standard output. Note that emissions are included as a source in  the Chemistry, and thus have the same timesteps. 


