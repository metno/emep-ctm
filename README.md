# EMEP MSC-W OPEN SOURCE CODE - 2017

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
`emep_Config_mod.f90` module. The aim is to restore `ModelConstants_ml.f90` to
its original role as a place to store some model-specific variables, and
to let `emep_Config_mod.f90` take over the operations of reading namelists
and setting the configuration variables. This work has only just started,
however, and so unfortunately the current release has configuration 
variables in both modules. This will be rectified in future release.

(The suffix `_mod` will also be used in future to represent module, rather
than ml as used prevously, since `_mod` seems to be gaining ground among
other fortran projects.)
