# EMEP MSC-W open-source release rv4.17, 2018


The EMEP MSC-W code released as rv4.15 (and with an emep-internal update
to rv4.16) was the result of many changes compared to previous versions,
and partly as a result of this we have found both bugs and areas of
improvement in the model. We summarise below the main changes made in
the rv4.17 (and previous interim rv4.16) versions.


# Land-cover issues

## Landuse_ml

 The rv4.15 code had several issues, especially when trying to run
 the  EMEP model in domains away from Europe (e.g. Asia). The new code
 improves initialisations, and can cope with the omission of the Euriopean
 landcover data.

 One issue remains: the code still expects two landcover files to be
 named in the config settings. If not using the default:

```
 LandCoverInputs%MapFile =
    'DataDir/Landuse/Landuse_PS_5km_LC.nc', 
    'DataDir/LandInputs_Feb2018/glc2000xCLMf18.nc',
```

One can use a dummy 2nd file:

```
 LandCoverInputs%MapFile =
    'DataDir/LandInputs_Feb2018/glc2000xCLMf18.nc',
    'DataDir/LandInputs_Feb2018/glc2000xCLMf18.nc',
```

## Global landcover

  NEW landcover file: glc2000xCLMf18.nc - set this in config_emep.nml.

  In the rv4.15 glc2000mCLM.nc file, deserts were included under the
  category 'BARE' instead of 'DE'. In principal this could have been
  fine, except that some 'BARE' land was also assigned where the GLC2000
  database had sparse vegetation. The EMEP model´s calculation of dust
  production was taking place only over the 'DE' regions, which were
  only found in the inner EECCA domain.

## Biogenics_ml

  The GetEuroBVOC routine in Biogenics_ml.f90 needed the landcover codes
   CF,DF,NF and BF to be defined. Those are defined by the European land-cover
   file Landuse_PS_5km_LC.nc. Users running for e.g. Asian domains use only
   the global land-cover file, which does not have these codes. The rv4.17 
   update fixes this issue.


# Radiation_ml

  As of version rv4.16 of the model, the radiation scheme of Weiss & Norman (1985)
  is used to calculate PAR levels. This gives a better treatment of diffuse and
  direct radiation compared to previous versions. 
  

# Chemistry - EmChem16a

EmChem16a is an update of the EmChem16 scheme presented in EMEP report 1/2017.

## Modifications

  Version rv4.15 of the EMEP model introduced a number of gas-aerosol
  reactions as part of the global study of Stadtler et al., 2018.
  Some of these gas-aerosol reactions have now been removed, either
  as they were negligable  (see Stadtler et al., 2018) or having a
  significant influence but causing model-performance degredation and
  being too uncertain to be reliable (NO2+aerosol). The main impact
  compared to rv4.16 is to raise NO2 concentrations to some extent.

  As part of the rv4.16 update, dry and wet deposition of N2O5, using
  the same rate as HNO3, was added to the chemical mechanism.

## bug fixes

  OH + HONO - rate coefficient had incorrect negative temperature dependence (small change)
  OD + H2O  - updated rate coefficient from IUPAC 2007

# config_emep.nml
  
    Pathes to all input files can now be steered by configuration settings.
    A short explanation for making own output fileds using 'USET' is given in the user guide.

    The n2o5Hydrolysi  method was incorrectly set in the 2017 open-source. It should be 'Smix',
    which is now the default in Config_module.f90.

# Timefactors

  For EECCA runs, the monthly time-factors from the LOTOS-EUROS model can now be
  activated, and seem to generate better results. 

  To use, set MonthlyNH3 = 'LOTOS' in config_emep.nml

# Tidy-ups

  A large number of tidy-ups have been done, and are underway. For example, the code
  used a lot of configutation variables (set via config_emep.nml) with the prefix
  USES\_. Each of these had to be named in the namelist settings. Now, the more 
  general user-defined variable USES covers many of these cases, with e.g.
  USES%FOREST_FIRES, USES%DEGREEDAY_FACTORS

# snow depth
  
  Modification of snow depth defintions: By default snow depth in the
  meteo fields is assumed to be in water equivalent, and multiplied by
  5 to give snow depth. For WRF metdata, the field named SNOWH is used
  (was SNOWNC in earlier versions); the WRF snowdepth is assumed to be
  directly in snow meters (i.e. not multiplied by 5).

# femis.dat

  "lonlat" type reductions has an additional country flag.

  The "lonlat" reductions can also be dectivated for individual emission
  files, using the configuration flag emis_inputlist(1)%use_lonlat_femis
  = F

# Radiation

  A height instead of level defintion is used. Should be more correct
  over mountainous terrain. Also in Southern hemisphere the wrong latitude
  was taken, now corrected.

# Acknowledgements

Some of the improvements made in the code stem from issues raised by EMEP users, either
via the github issues sides or personal communication. Thanks are due to 
John Johansson and Robert Bergström (Chalmers, Sweden), Roy Wichink Kruit (RIVM),
Paul Hamer (NILU), and Massimo Vieno (CEH, Scotland) for such feedback. 


# References

Stadtler, S., Simpson, D., Schröder, S., Taraborrelli, D., Bott, A., and Schultz, M.: Ozone Impacts of Gas-Aerosol Uptake in Global Chemistry Transport Models, Atmos. Chem. Phys., in press, 2018

Weiss, A. and Norman, J. M.: Partitioning Solar-radiation into Direct and Diffuse, Visible and Near-infrared Components, Agricultural and Forest Meteorology, 34, 205–213, doi:10.1016/0168-1923(85)90020-6, 1985.


# EMEP MSC-W open-source release rv4.15 - 2017


The source code in this 2017 release is essentially that used to model
air pollution in the year 2015 (rv4.15), as reported in EMEP MSC-W Status
Report 1/2017 (www.emep.int).  For basic documentation of the model see
the file CITATION.txt which is released with the code.  A major difference
in the model setup this year is a finer resolution, namely 0.1x0.1
degree in the horizontal and 34 layers in the vertical.  Furthermore,
this model version has had several major coding and data-input changes
since the 2016 release, and in many ways should be seen as an interim
version of the model. Here we just list some specific issues associated
with the code and model usage.


## 34-layers

The use of 34 layers rather than 20 is very new to the EMEP model,
and some aspects of the model formulation may need to be revised if we
continue with this setup.  In particular the dry-deposition formulation
was originally designed for thicker (ca. 90 m) lower-layers, and may
need modification to cope with the new structure.

The revised model still works however with the 20 layer meteorology as
provided with previous releases, and indeed is still flexible in terms
of allowing additional layers within the 20-layer structure (through
interpolation) or for use with e.g. WRF meteorology.


## Shipping emissions

The new IMO (International Maritime Organization) regulation on sulfur
emissions from international shipping, which came into effect in January
2015, has led to a significant reduction in sulfur emissions within
the so-called Sulfur Emission Control Areas (SECAs) in Europe, i.e. the
Baltic Sea and the North Sea.  The MACC-TNO-III data set for 2011, which
was used for EMEP reporting until last year could thus not be used any
longer, 2011 SOx emissions are significantly larger than those of 2015.

A new data set (see Status Report 1/2017) has been made available to us
for 2015 by the Finnish Meteorological Institute (FMI). It is based on
accurate ship position data from AIS (Automated Identification System)
and also takes into account the new IMO regulations. As part of the
Copernicus Atmospheric Monitoring Service (CAMS), ship emissions will
be calculated by FMI also for years after 2015, and it is assumed that
data for 2016 will arrive in time for next year's EMEP reporting.

Unfortunately, this situation results in an inconsistency in data-sources
between the MACC-TNO-III data and the FMI data. Users should be aware of
this and make appropriate choices when running years other than 2015.

Implementation of these new emissions into the EMEP model was rather
hurried, with the new emissions in a special format, and with specific
code changes to cope with this format. For future releases of the model
we will re-format the emissions to comply with the usual EMEP formats,
and eliminate special coding as far as possible.

## Code structure

Some code from `ModelConstants_ml.f90` has been moved into a new
`Config_module.f90` module. The aim is to restore `ModelConstants_ml.f90` to
its original role as a place to store some model-specific variables, and
to let `Config_module.f90` take over the operations of reading namelists
and setting the configuration variables. This work has only just started,
however, and so unfortunately the current release has configuration 
variables in both modules. This will be rectified in future release.

