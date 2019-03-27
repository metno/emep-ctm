.. _`ch-submitarun`:

Setting the input parameters
============================

In this chapter we provide detailed information on how to set the parameters of the regional EMEP/MSC-W model.

In general the parameters and pathes to the input files are all set in the configuration file ``config_emep.nml`` (a fortran namelist).

``config_emep.nml``
-------------------

The model has a namelist system. It is possible to set different
constants and flags for running the model. The constants and flags
themselves are defined in ``Config_module.f90``, while they are set in the
namelist file under ``ModelConstants_config`` parameter. Some of these
are briefly explained in :numref:`ch-inputfiles`. Model gets information
about running for special cases from this file. The datasets provided
are for the EMEP grid EECCA.

The different parameters for the model run are set in the
``config_emep.nml`` file. In the very beginning of this, the section
``INPUT_PARA`` has all these variables including the link to the
meteorology data. The trendyear can be set to change the boundary
emissions for earlier and future years, see the modules
``BoundaryConditions_ml.f90`` to understand
better what the trendyear setting does. The default setting is the
meteorological year you are running for, in this case 2015. The
runlabel1 option sets the name of the different output NetCDF files, see
:numref:`ch-output`. The ``startdate`` and ``enddate`` parameters are set for the
time period you want the model to run (YYYY,MM,DD), and you need
meteorology data for the period, as shown below:

.. code-block:: text
  :name: config-emep
  :caption: Basic namelist example; ``config_emep.nml`` extract.

  &INPUT_PARA
    GRID = 'EECCA',
    iyr_trend = 2015,
    runlabel1 = 'Base',
    runlabel2 = 'Opensource_Setup_2018',
    startdate = 2015,01,01,00,
    enddate = 2015,12,31,24,
  &end
  &Machine_config
   DataPath(1) = '../input', ! define 'DataDir' keyword
  &end
  &ModelConstants_config
   meteo = '../meteoYYYY/GRID/meteoYYYYMMDD.nc',
   DegreeDayFactorsFile = 'MetDir/DegreeDayFactors.nc',
  !------------------------------
   EmisDir = 'DataDir/GRID',
   emis_inputlist(1)%name= 'EmisDir/gridPOLL', !example of ASCII type
  !--------Sub domain x0, x1, y0, y1
    RUNDOMAIN = 36, 100, 50, 150, ! EECCA sub-domain
  &end

In the extract above, the model is run for the period 1 January to 10 January
2014 and the trend year used is 2014. Output files will be stored with
the name 'Base' and the meteorological correspond to the 'EECCA' grid.

It is possible to run the model on a smaller domain than the full
regional model domain, as defined by indexes :math:`x` and :math:`y`.
For the 'EECCA' grid  :math:`x=1,\ldots,132; y=1,\ldots,159`\ .
To set a smaller domain, use ``RUNDOMAIN`` variable in the ``ModelConstants_config``
namelist to indicate the sub-domain indexes. In the config_emep extract above,
``RUNDOMAIN`` defines a subdomain with :math:`x=36,\ldots,100; y=50,\ldots,150`\ .

Using and combining gridded emissions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The gridded emission files are controlled via the ``config_emep.nml``
file. Each file is assigned as one set of values for ``emis_inputlist``.
:numref:`emis-config` line 1 includes an ASCII emission file, where the
keyword ``POLL`` will be replaced by the model by all the
emitted pollutants (according to the names defined in ``CM_EmisFiles.inc``).
An additional NetCDF emission is included in line 2.

Now all emissions from both ASCII file and NetCDF file will be used. In
practice some countries might be counted twice. Therefore some new data
can be included in the ``emis_inputlist``, to specify which countries to
keep or to avoid. :numref:`emis-config` lines 3--4
will include only 'NO', 'SE' and 'FI' from the first file (ASCII), and
take all countries except 'NO', 'SE' and 'FI' from the second file
(NetCDF).

Sets of countries can in principle be defined; for now only the set
'EUMACC2' is defined.

.. _`emis-config`:

.. code-block:: Fortran
    :caption: Mixed emission configuration example.
    :linenos:

    emis_inputlist(1)%name = '/MyPathToEmissions/emislist.POLL',
    emis_inputlist(2)%name = '/MyPathToEmissions/Emis_GLOB_05.nc',
    emis_inputlist(1)%incl(1:) = 'NO','SE','FI',
    emis_inputlist(2)%excl(1:) = 'NO','SE','FI',
    emis_inputlist(1)%PollName(1:2) = 'voc','sox',


It is also possible to restrict the number of pollutants from each of the files.
If not all pollutants from ``CM_EmisFiles.inc`` are to be read, one can specify a list of pollutants to be included
using "PollName". For instance in the example above , the last line specifies that emissions will include only VOC 
and SOx from the file defined by emis_inputlist(1)%name. 
If PollName is not specified at all, all pollutants are included (therefore all pollutants 
from emis_inputlist(2)%name will be included).
The specified pollutants must already be defined in ``CM_EmisFiles.inc``.
It is possible to disregard the "lonlat" reductions introduced by ``femis.dat`` for specific emissions. To do this use the "use_lonlat_femis" flag.
Example: switch off emissions covering one region from ``Emis_GLOB_05.nc`` as specified by femis, and replace the emissions in that data using ``emislist.POLL``

.. code-block:: Fortran
    :caption: Do not take into account the lines starting with lonlat in femis.dat for ``emis_inputlist(1)%name``.
    :linenos:

    emis_inputlist(1)%name = '/MyPathToEmissions/emislist.POLL',
    emis_inputlist(1)%use_lonlat_femis = F,
    emis_inputlist(2)%name = '/MyPathToEmissions/Emis_GLOB_05.nc',


An alternative way of combining overlapping emissions, is to use a "mask" approach. This is typically used, when the emissions of one city is known in more details, but one wish to include default regional or global emissions elsewhere.
The city emissions are used to set the mask, and that mask is then in turn used by the second emission sources to turn off emissions within the city.

.. code-block:: Fortran
    :caption: Mixed emission using mask configuration example.
    :linenos:

    emis_inputlist(1)%name = '/MyPathToEmissions/emis_local.POLL',
    emis_inputlist(1)%set_mask = T,
    emis_inputlist(2)%name = '/MyPathToEmissions/Emis_GLOB_05.nc',
    emis_inputlist(2)%use_mask = T,

There are some details one should take into account: the order of the "set_mask" and "use_mask" matter. The mask should be set before it is used; "before" meaning that the indice of "emis_inputlist" is lower. Also yearly emissions are treated before monthly emissions, therefore one should not use monthly emissions to mask yearly emissions.
The mask is set for a given position if emissions at that point are larger than 1.0e-20. If emissions are zero at some point the mask will not be set for that point, and regional emissions will be included there.
There is only one mask. Several emissions files can set and use the mask.

config: Europe or Global?
~~~~~~~~~~~~~~~~~~~~~~~~~


The EMEP model has traditionally been run on the EMEP grid covering
Europe, and using meteorology from the ECMWF IFS model. In this environment,
we typically set several configuration variables to make use of
Euro-specific data. In other regions it is better to make use of the
model's 'global' settings, which will ensure better handling of
vegetation (eg LAI) changes, convection, and various emission
settings.

.. code-block:: text
  :caption: Typical European/EMEP settings.

  USES%DEGREEDAY_FACTORS = T,    ! though F is okay too
  USES%PFT_MAPS = F,
  USES%MonthlyNH3 = 'LOTOS', ! Better monthly profile, for Europe only!
  USES%CONVECTION = F, 
  USES%EURO_SOILNOX = T, ! diff for global + Euro

.. code-block:: text
  :caption: Typical non-European/global settings.

  USES%DEGREEDAY_FACTORS = F
  USES%PFT_MAPS = T,    ! PFT LAI tests
  USES%MonthlyNH3 = '-', ! Better monthly profile, for Europe only!
  USES%CONVECTION = T
  USES%EURO_SOILNOX = F, ! diff for global + Euro runs
  USES%GLOBAL_SOILNOX = T, ! diff for global + Euro runs

The ``DEGREEDAY_FACTORS`` setting triggers the use of degree-days in
controlling SNAP2 emissions. This requires pre-processed files of heating degree days.
Such files (DegreeDayFactors.nc) can be produced from any meteorology, but the difference in
results even for Europe is not too significant. In other regions of the
world emissions from SNAP2 may not be as dependent on degree-days as in
Europe, and so this setting should probably be false.

``PFT_MAPS=T`` triggers the use of a global file which provides monthly
variations in leaf area index (LAI) for different vegetation types. This
controls deposition and biogenic VOC emission parameters. In Europe, a
simpler latitude-dependent system is used (based upon DO3SE), and so
``PFT_MAPS`` should be set F.

``MonthlyNH3='LOTOS'`` is also only relevant for European simulations; and
indeed any non-European runs are better off with monthly emissions for
that particular area.

``CONVECTION`` is difficult. In principle, all models runs should use the T
setting, but for Europe we find it degrades the model results too much
and we use F. The problem is likely that the sub-grid processes behind
convection are so complex and the paramererisation is very uncertain.
Note also that in Config_module we have the default setting ``CONVECTION_FACTOR=0.33``,
which may be changed to allow more or less influence of this variable.



Input file
----------

Vertical_levels.txt:
The number in Vertical_levels.txt correspond to the "A" and "B" coefficients of the hybrid (eta) coordinates (P=A+B*Psurf).

Close to the surface, A should be small, and higher up we should use pressure levels. Then there is a gradual transition from surface to pressure levels.

For the case of the Vertical_levels20.txt; the levels are identical to the sigma layers we used originally. Sigma levels are a special case of hybrid levels.



Base run
--------

This is an example of a minimum ``modrun.sh`` script to run the model.

.. literalinclude:: modrun.sh
    :caption: Minimum ``modrun.sh`` example.
    :language: bash

This bash shell script is designed so that users can easily adapt it to
fit their needs. It contain the minimum information required to run the
EMEP/MSC-W model. The script should be self-explanatory. The paths 
for meteorology and all other input data are defined on the model configuration file,
``config_emep.nml``. You need to set the right paths for the input directories.

It is possible to run the model on a smaller domain than the full
regional model domain, as defined by indexes :math:`x` and :math:`y`.
For the 'EECCA' grid  :math:`x=1,\ldots,132; y=1,\ldots,159`\ .
To set a smaller domain, use ``RUNDOMAIN`` variable in the ``ModelConstants_config``
namelist to idicate the sub-domain indexes. In :numref:`config-emep`,
``RUNDOMAIN`` defines a subdomain with :math:`x=36,\ldots,100; y=50,\ldots,150`\ .

To run the model, the correct path to the EMEP/MSC-W model code has to
be set (``mpiexec path_to_the_modelcode/Unimod``).

It is recommended to submit the script as a batch job. Please check the
submission routine on the computer system you are running on.
In the newer model versions (since 4.0) the number of nodes is set
automatically from what is asked for when submitting a job.
The approximate time and CPU usage is described in :numref:`sec-compinf`.

When the job is no longer running or in the queue, it is either finished
or has crashed for some reason. If the model run crashed, an error
message will give information on what was missing or wrong in the
routine. If the run was successful, a message

.. code-block:: text

    ++++++++++++++++++++++++++++++++++++++++++++++++
    programme is finished

will be stated at the end of the log file, before printing the
``Timing.out`` file.

The model results will be written to this same
directory. Please make sure there is enough disk place for the model
results. For more information about the model result output files,
see :numref:`ch-output`.

If for some reason the model crashed, please check both the log and the
error file for any clue of the crash. After fixing the problem the job
can be submitted again. If the model has crashed, then the links to the
input data are not removed.

The script can also be submitted interactively, and either have the
output written to the screen or to named error and output log files. The
variables wanted in the output are specified in the
``OutputConcs_config``, ``OutputDep_config`` and in the ``OutputMisc_config``
parameters respectively for surface concentrations, depositions and some
miscellaneous outputs.

Separate hourly outputs
-----------------------

The ``Base_hour.nc`` and ``Base_uEMEP_hour.nc`` files can become very lkarge. It is possible to split them into one file per day by adding 
the keyword ```HOURLYFILE_ending    = 'JJJ.nc' ```, in the configuration file. ```JJJ```, will be automatically replaced by the corresponding day of the year (i.e. a number from 1 to 366 giving for instance
``Base_hour_001.nc``). The full date as in ``Base_hour_20180101.nc`` can be obtained by defining ```HOURLYFILE_ending    = 'YYYYMMDD.nc' ,```. 

.. _`sec-nesting`:

Nesting
-------

The boundary conditions needed for EMEP MSC-W model is provided with the
input data. The model can read Boundary conditions data from other
models as well. These data has to be in NetCDF format. The boundary
conditions needed for EMEP MSC-W model is provided with the input data.
The model can read Boundary conditions data from other models as well.
These data has to be in NetCDF format.

Different Nesting modes are:

*  read the external BC data only,
*  produce EMEP BC data from the simulation,
*  read the external BC data and produce EMEP BC data,
*  using the default EMEP BC data from the input data directory and
   write out EMEP BC at the end of the simulation,
*  read the external BC data only in the beginning of the simulation,
*  read external BC at the beginning of the simulation and write out
   EMEP BC at the end of the simulation.

These options are controlled by the ``MODE_READ`` and ``MODE_SAVE`` variables
in ``Nest_config`` namelist, in ``config_emep.nml`` file. The mode options are:

``MODE_READ``
    'NONE'
        do nothing (default).
    'START'
        read at the start of run.
    'FORECAST'
        read at the start of run, if the files are found.
    ‘NHOUR’
        read at given ``NHOURREAD`` hourly intervals.
        ``NHOURREAD`` is set in ``Nest_config`` and should be an integer fraction of 24.
    ‘MONTH’
        read at start of each month.
``MODE_SAVE``
    'NONE'
        do nothing (default).
    'END'
        write at end of run.
    'FORECAST'
        write every ``OUTDATE(1:FORECAST_NDUMP)``.
        ``OUTDATE`` and ``FORECAST_NDUMP`` are set in ``Nest_config``.
    'NHOUR'
        write at given `NHOURSAVE` hourly intervals.
        ``NHOURSAVE`` is set in ``Nest_config`` and should be an integer fraction of 24.
    ‘MONTH’
        write after each month (used for checkpoint/restart for instance).

BC data are read from the file defined by ``template_read_BC`` in ``Nest_config``. 
In addition, the IC data (entire 3D domain) will be set at start of run from the file defined by ``template_read_3D`` in ``Nest_config`` (except if  ``MODE_READ`` = 'NONE').

Write BCs from EMEP MSC-W model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:numref:`nest-write-config` shows an example to write every 3 hours into
daily Nest/BC files. Output file name is defined by ``template_write`` ('BC_YYYYMMDD.nc'),
where ‘YYYYMMDD’ is replaced by corresponding date.

All advected model variables will be written out for the sub-domain defined
by ``out_DOMAIN`` (\ :math:`x=60,\ldots,107; y=11,\ldots,58`\ ).
If no ``out_DOMAIN`` is given, the entire model rundomain will be written out.

.. code-block:: text
    :name: nest-write-config
    :caption: Write BCs configuration example.

    &Nest_config
      MODE_READ        = 'NONE',          ! do not read external BC
      MODE_SAVE        = 'NHOUR',         ! write BCs
      NHOURSAVE        = 3,               !  every 3 hours
      template_write   = 'BC_YYYYMMDD.nc' !  to your (daily) BC output file
                                          !  (8 records per file)
      !-------- Sub domain for write modes
      out_DOMAIN       = 60,107,11,58,    ! istart,iend,jstart,jend
    &end
    
Reduce the size of BC files
___________________________

The size of the files obtained in a nesting configuration can be very large if the out_DOMAIN is large.
If the inner domain is known in advance, only the part matching exactly the part needed to construct the BC can be stored.
Define ``MET_inner`` in ``&Nest_config``, which should be a link to any metdata of the inner grid;
  it will only be used to define the projection parameters of the inner grid (i.e. dates and other content do not matter).

.. code-block:: text
    :name: nest-write-inner
    :caption: Inner domain options for nested BC output example.

    &Nest_config
    [...]
      MET_inner = 'inner_domain/wrfout_d03_YYYY-MM-DD_00:00:00'
    &end
        
Note that the file will have the same dimensions, but zeros are put into the unused parts.
The NetCDF internal compression will take care of reducing the actual size, as measured by used disc space.

If a BC file has been created using the ``MET_inner`` method, it cannot be used for initializating concentrations at the start of the run. A separate file has to been created for initializations. This file can then be used by the inner grid by defining ``template_read_3D`` in ``config_emep.nml``.


Read BCs from EMEP MSC-W model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:numref:`nest-read-config` shows an example to read every 3 hours from
the Nest/BC files created previously by running :numref:`nest-write-config`.
Please note that the model sub-domain for a nested run is set by ``RUNDOMAIN``,
as shown in :numref:`config-emep`.

.. code-block:: text
    :name: nest-read-config
    :caption: Read BCs configuration example.

    &Nest_config
      MODE_READ        = 'NHOUR',         ! read external BC
      NHOURREAD        = 3,               !   every 3 hours
      template_read_BC = 'BC_YYYYMMDD.nc' !   your (daily) BC input file
      MODE_SAVE        = 'NONE',          ! do not write BCs
    &end

Read external BCs
~~~~~~~~~~~~~~~~~

Reading BCs from a different model is more involved than the previous example.
The vertical axis and variables in the file need to be mapped to the
corresponding model variables.

:numref:`nest-mybc-config` shows an example to read every 3 hours from an external
BC file. The model will read 3 variables from ``MyBC.nc``: |O3|, NO, and |NO2|.
The maping between the ``MyBC.nc`` variables and the corresponding model variables
is defined in the ``ExternalBICs_bc`` namelist.

.. code-block:: text
    :name: nest-mybc-config
    :caption: External BCs configuration example.

    &Nest_config
      MODE_READ        = 'NHOUR',         ! read external BC
      NHOURREAD        = 3,               !   every 3 hours
      template_read_BC = 'MyBC.nc'        !   from your BC input file
      MODE_SAVE        = 'NONE',          ! do not write BCs
    &end
    &ExternalBICs_config
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

Unit conversions are delegated to the ``Units_ml.f90`` module.
The supported units are: |ugm3|\ , |ugSm3|\ , |ugNm3|\ , |ugCm3|\ , ppb,
mixing ratio (mol/mol) and mass mixing ratio (kg/kg).

If the BC data has different units,
either convert them into one of the above mentioned units in pre-processing
or add the respective conversion factor in the module ``Units_ml.f90``.


Source Receptor (SR) Runs
-------------------------

The EMEP/MSC-W model can be used to test the impact of reduced emission
of one or more pollutants from a particular country or a number of
countries. Such runs are called "Scenario runs". They are used for source-receptor calculations.

Emission factors for reduced emissions of pollutants from different
sectors and countries can be defined in the input file called
``femis.dat``, which can be found in the downloaded input data directory,
see section :numref:`sec-femis`.

.. code-block:: text
    :name: base-femis
    :caption: ``femis.dat`` for a base run.

    Name  7  sox  nox  co   voc  nh3  pm25   pmco
    27    0  1.0  1.0  1.0  1.0  1.0  1.0    1.0

This base run example means that there are (1.0), no emission reductions
of |SOx|\ , |NOx|\ , CO, VOC, |NH3|, |PM25| and |PMco| from all sectors in the UK.

-  The first column of the second line represents the country code. (27
   is the code for UK.) The codes for all countries can be found in
   Fortran module ``Country_ml.f90`` or at (http://www.emep.int/grid/country_numbers.txt). 
   The country code must be the same as in the emission files for the given country. Some
   countries and areas are divided into sub-areas in the emission files.
   In this case, one line for each sub-area has to be included into the
   ``femis.dat`` file. Countries and areas where emissions are given for
   sub-areas include the Russian Federation, Germany and all sea areas.
   "0" means all countries.

-  The second column of the second line represents the sector and "0"
   means all sectors. Here one can write the appropriate sector code if
   the emission is reduced only from a specific sector. The description
   of each sector can also be found in the Fortran module ``EmisDef_ml.f90``.

-  The columns under the pollutant names show the emission factors for
   the given pollutants. For example, 0.7 would mean 70% of the original
   emission, thus 30% reduction.

-  The number ("7") following the first text ("Name") in the first line
   gives the number of pollutants treated in the file.

An example of ``femis.dat`` file describing 50% reduced emission of |NH3|
from sector 10 (the emission from agriculture) in the UK is shown in
:numref:`reduction-femis`.

.. code-block:: text
    :name: reduction-femis
    :caption: ``femis.dat`` for 50% |NH3| reduction from sector 10 over UK.

    Name  7  sox  nox  co   voc  nh3  pm25   pmco
    27   10  1.0  1.0  1.0  1.0  0.5  1.0    1.0


Instead of entire countries, reductions can also be specified by coordinates (and combined with country reductions).
The line with coordinate corrections must start with the keyword ``lonlat``. The coordinates are given in longitude latitude (min and max and the coordinates of the centre of the gridcells are tested. Gridcells are either entirely included or entirely reduced, never cut into smaller parts).


.. code-block:: Fortran
    :name: femis
    :caption: ``femis.dat`` example.
    :linenos:

    Name                          7  sox  nox  co   voc  nh3  pm25  pmco
    17                            0  1.0  1.0  1.0  1.0  1.0  0.5   0.5
    lonlat 3.3 7.2 50.7 53.5   17 0  1.0  1.0  1.0  1.0  0.0  1.0   1.0


In :numref:`femis`, country with code 17 (NL) will reduce |PM25| and |PM10| emissions by half for all sectors.
Emissions of |NH3| from country with code 17 only, will be removed from the rectangle with longitudes between
3.3 and 7.2 degrees East, and between 50.7 and 53.5 degrees North. Use zero (0) as country code to specify that emissions from all countries should be reduced.


New emission format
-------------------
A new more general and (hopefully) easy to use format for emissions has been introduced. It is still in a developing phase, so changes and errors may occur.

In the new format, emissions are organised in a number of files (Emis_sourceFiles(i_file)), each files containing a number of sources (Emis_sourceFiles(i_file)%source(j_source)).
For now the main constraint is that a source is any 2D field (possibly+time).
The file must have a ‘lon’ and a ‘lat’ variable, showing longitude and latitudes of each grid point. ‘lon’ and ‘lat’ must be 1D variables if the projection is ‘lon lat’, 2D otherwise.

The file and sources can be characterized by a set of variables. In genereal these variable can be set by and in order of increasing priority:
  1. Default value
  2. Global attribute read in the netcdf file
  3. Variable attribute read in the netcdf file
  4. Value set for Emis_sourceFiles(i)%XXX in config_emep.nml
  5. Value set for Emis_sourceFiles(i)%source(s)%XXX in config_emep.nml
  
Exception to the priority rule are:
  * maskID cannot be set by attributes in the netcdf file
  * the file and source factor are on top of each other, not replaced
  * boolean parameters (like apply_femis), are used as "and" (i.e. if any is false, the result is false) 

List of file attributes (default in parenthesis):
  - filename (‘NOTSET’) Name of the file (with path)
  - projection (‘lon lat’) Only three categories ‘lon lat’, 'native' or any other (for example ‘Lambert ‘or ‘Stereographic’ would give the same result). 'native' means that emissions are given in the same grid as the model grid and the data is not interpolated (use it if you can).
  - grid_resolution (an approximate value is computed from the lon and lat, if no value is given) It does not need to be exact (cannot be exact on a sphere anyway!). This grid_resolution steers the interpolation algorithm; A large value will force the code to subdivide each emission gridcell in large number of pieces, that are assigned to the model grid. Larger values means smoother interpolation, but more cpu time. 
  - periodicity (‘time’) How often the values are updated. Can be ‘yearly’, ‘Monthly’, ‘hourly’ or ‘time’. ‘hourly’ or ‘time’ means that the time as defined in the netcdf is used to define when to fetch a new record. The timestamp must correspond to the end of the time period of validity. For ‘yearly’ monthly timefactors are applied, if a sector is defined. For ‘monthly’ and ‘yearly’, an hourly timefactor is applied if a sector is defined. For ‘hourly’ or ‘time’, no additional timefactors are applied. 
  - factor (1.0) multiplicative factor for all sources in the file
  - apply_femis (true) whether to apply the femis reductions to the sources of this file.
  - include_in_local_fractions (true) whether to take sources from this file into account for the local fraction calculations
  - units ('NOTSET') will be used as default for sources units if set.
  - country_ISO ('NOTSET') will be used as default for sources units if set.
  - sector (-1) will be used as default for sources sector if set.
  - species ('NOTSET') will be used as default for sources species if set.
  - mask_ID ('NOTSET') the name of the mask, if you want to apply one.
  - mask_ID_reverse ('NOTSET') the name of the mask, if you want to apply one in the complementary region.

List of source attributes:
  - varname (‘NOTSET’) The name as used in the netcdf file
  - species (‘NOTSET’) Either one of the emission group species, as defined in CM_EmisFile.inc (generally sox, nox, pm25, pmco, nh3, co, voc) 
  - factor (1.0) multiplicative factor. Can be used to change units to model definitions.
  - units (‘mg/m2/h’) Units *after* the factor multiplication. Comes on top of the file multiplicative factors and possibly other factors.
  - country_ISO (‘N/A’) the country code, as defined in Country_mod.f90 (for example ‘FR’ for France). ‘N/A’ is a valid code, but it does not correspond to any country.
  - apply_femis (true) whether to apply the femis reductions to this source.
  - include_in_local_fractions (true) whether to take this source into account for the local fraction calculations
  - mask_ID ('NOTSET') the name of the mask, if you want to apply one.
  - mask_ID_reverse ('NOTSET') the name of the mask, if you want to apply one in the complementary region.


The idea is that only variables that clearly are required in a specific context need to be set; if the value can be inferred from other information, the code should do it.
Depending of the type of source, not all variables are used.

Note about species: These can be interpreted in one of three categories
  1. emitted species (nox,sox,pm25 ...) with sector (1...11) (“sector species”)
  2. individual species (SO2, NO, NO2, ...) with sector. The species MUST be one of the splitted species. These will be treated as one of the “sector species”  from 1. (but not splitted of course). Careful with units, it follows the same rules as “sector species”; molecular weight for SO4 for example is considered “as SO2”.
  3. individual species (SO2, APINENE, O3 ...) without sector (<=0, or not specified). No timefactors, vertical realease heights or splits are applied.
   In this case the emissions are summed up in setup_rcemis (not in EmisSet)
   
Masks
-----

Typically you have got fine scale emissions for a small region of interest, a city for instance. You may want to remove that area from the coarse scale emissions, and replace it with your own. The mask allows you to define a specific region (the mask).
To define which gridcells to include in your local region, you must find a suitable variable that shows the region of interest. It could be for example the PM emissions in your local area. 

A "mask" can be defined for instance with:

.. code-block:: Fortran
    :caption: Define a mask example
    :linenos:

    EmisMask(1)%filename = '/mypath/myfile.nc' !name of the netcdf file to read from
    EmisMask(1)%cdfname  = 'London_PM'  !name of the variable to read from the file
    EmisMask(1)%ID       = 'LONDON'  !the name you give to that mask
    EmisMask(1)%threshold = 1.0E-10 !the mask is set at any point larger than the threshold

    
Several masks can be defined. Each mask is identified by their "ID". If you want to include in the region also the gridcell whicgh are zero, you can set the threshold slightly negative (-1.0E-10), to include the entire region covered by the variable (otherwise zero values would be defined equivalently to outside of region).

A mask defines only a region. It is not directly related to any pollutant. 

The masks defined here, will also be applied on files from emis_inputlist (old format), if use_mask is set. It is however not possible to set masks by both systems simultaneously. In the old format only one mask can be used at a time. It will be the reunion of all masks produced by the system above (the ID is meaningless and cannot be specified in old format).

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
See also in ``Par_ml.f90`` for details.

The main timestep parameter ``dt_advec`` can be set manually in config_emep.nml (in seconds, real number).
It must be a factor of 3600.

``JUMPOVER29FEB`` : if set to T , will not treat the 29th of February even for leap years.
(Useful if that date is missing, for instance in Climate runs).

``NETCDF_DEFLATE_LEVEL``: The level of compression used in the NetCDF output files (integer).
Negative values means netcdf3 format.

Defining own sectors
--------------------

Emissions can be assigned to a sector. A sector defines three properties:

1. Emission height release distribution
2. Split into species
3. Timefactors

The set of Emission heights available is defined in the file "EmisHeights.txt".  The set of splits in files "emissplit.defaults.POLL" and "emissplit.specials.POLL". The timefactors are defined in "MonthlyFac.POLL" "DailyFac.POLL" and "HourlyFacc.INERIS".

Which height/split/timefac is chosen for a given sector is defined through a mapping system. There are two predefined mapping you can choose from: SNAP and GNFR. The mapping are simply three tables (one dimensional array). Those mappings are defined in "EmisDef_mod.f90" in the arrays "XXX_sec2hfac_map", "XXX_sec2sfac_map","XXX_sec2tfac_map", where "XXX" is the name of the mapping (SNAP or GNFR).
For example the mapping "GNFR_sec2hfac_map = (/1,3,2,4,6,7,8,8,8,9,10,10,5/)", means that in the GNFR convention, the sector 13 is mapped to the fifth emission height in "EmisHeights.txt".
You can add both more emission heghts in "EmisHeights.txt" and access them by changing the maps. (The maps cannot be set by config_emep.nml for now because it is not guaranteed that everything will work when you change the wrong things).
You can define a new mapping for example using the "TEST" mapping (also in EmisDef.f90). To switch to "TEST" for those mappings set "USE_SECTOR_NAME='TEST'" in "config_emep.nml".


