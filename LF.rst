
Local Fractions
===============

The Local Fraction method allows to give information about where the pollutants come from. There are two main "branches": tracking of primary particles, and "generalized LF", which also includes chemical transformations between species. The values for the local fractions will be outputted in separate files (file names including  "_LF_"). 




Local Fractions for primary particles
-------------------------------------

Primary particles are simplest; we can imagine that pollutants from each source are tracked separately (it is what we do in the code). What is included as a source is defined by the settings in config_emep.nml.

Mainly two types of sources (you can have both type together in the same run):

* grid to grid (aka relative), which considers each individual gridcell separately. To keep the amount of data reasonable, the pollutants are tracked only up to a given distance (10 gridcells in each direction for example). ``lf_src(1)%type='relative'``, default.
* countries, possibly limited to a specific sector. (SR type runs). Specify ``lf_src(1)%type='country'``


Below an example of grid to grid LF config_emep.nml settings, which is the type used by the urban EMEP, uEMEP.


.. code-block:: Fortran
    :caption: Local Fractions flag example

    USES%LocalFractions = T, ! T for computing Local Fractions, F otherwise
    !Local Fractions frequency of output (separate file for each). Can be any of: YEAR, MONTH, DAY, HOUR, HOUR_INST 
    lf_set%YEAR = T, !average value for full run in output
    lf_set%Nvert = 14, !How many vertical level to include in treatment. Should be higher than highest emissions
    lf_set%dist = 5,  !how far the neighbors can be in each direction (NB: high cost for large dist)
    
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

The local fractions values given are for the lowest level; they are not interpolated to 3 meter using the "cfac" factor. 

Tip: if you want to look at the grid-to-grid fields in the netCDF file with ncview; press first on any regular 2- or 3-dimensional field. (otherwise ncview may put your zoom to X100 or so!). Then choose axes lon-lat, then press x_dist or y_dist to set them to zero or small.

"LF_GUI_LF.py" https://github.com/metno/emep-ctm/tree/tools can be used to visualize the grid-to-grid fields.



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

The results for country types outputs are interpolated to 3 meter height, using the "cfac" factor.
    
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
    lf_spec_out(10)%name='EUAOT40_Crops',
    lf_spec_out(11)%name='pm25',
    lf_spec_out(12)%name='PM_WATER',

    lf_set%MDA8 = T, !special: make AvgMDA8_6month and SOMO35 NB: requires that O3 is outputted too
    !Note that SOMO35 has a different defintion than in the "base" file output: Mainly in the "base" file SOMO35 is computed 
    !excluding the 8 hours periods that overlap with midnight. Also the "base" output include the 8 hours periods every 20 
    !minutes (independently of dt_advec!), while every hour is used for the LF output. 
    !The LF definition of SOMO35 follows similar rules as the more official and detailed MDA8. 


Miscellaneous
-------------

Even if we consider pm25 as "primary", there is some "chemistry" going on: the different pm species are "new" or "age" and will be transformed gradually from new to age. This is taken into account also in the LF for primary particles (the two classes are tracked separately).

For POD and AOT outputs, the name must match one of the output names used in the regular output (not local fraction output). 

'BVOC', 'DMS', 'BC', 'STRATOS', 'INIT' are special "countries".

'ALL' is a predefined group, including all countries.

When using nesting for restart, the local fractions will also be stored if activated. This will only work for "country" style LF, and only if the grid is identical before and after restart.

In order to compare with a brute force method, you can reduce a country emission with 1%. For better agreement use (in both runs): ``ZERO_ORDER_ADVEC = T`` , ``AERO%EQUILIB='MARS'``, ``AERO%EQUILIB_WATER='MARS'``
