# Open Source EMEP/MSC-W model

## Welcome
The EMEP models have been instrumental to the development of
air quality policies in Europe since the late 1970s,
mainly through their support to the strategy work under the
[Convention on Long-range Transboundary Air Pollution][CLRTAP].
In the 1990s the EMEP models became also the reference tools for
atmospheric dispersion calculations as input to the Integrated Assessment Modelling,
which supports the development of air quality polices in the European Union.

Through this web page the latest EMEP/MSC-W model version used for the
[EMEP status reporting of the year 2016][publ2016] - rv4.10 - is released,
together with a set of input data,
an updated user guide and a full year model results for the year 2014,
under [GPL license v3][GPLv3].

You are recommended to read all chapters of the
[EMEP/MSC-W model User Guide][guide]
before you start downloading data from the following link:<br/>
[OpenSource rv4.10 (201609)][rel410].

This contains the following set of information:

* a complete set of 'input data' to allow for model runs for year 2014
* the open source 'model code' of the EMEP/MSC-W model version rv4_10
* 'model results' for the year 2014 for comparison of a successful run

The standard input files for the EMEP model are provided in two different formats,
netCDF or ASCII.
We presently follow as much as possible the [netCDF CF conventions][netCDF_CF]
for both input and output data.
More details about the input files are described in Chapter 2 of the
[EMEP/MSC-W model User Guide][guide].

The tar file of the 'Open source model code' includes the model source code,
a run script, a makefile and a copy of the gpl license.
More details about the model code is given in the User Guide.  

The main output files are in netCDF ([netCDF CF conventions][netCDF_CF]) format
and some additional output files are in ASCII format.
Details about the output files are described in Chapter 3 of the User Guide.

How to submit a run is described in detail in the chapter 4,
'Submitting a Run' of the User Guide.

The EMEP/MSC-W model is very flexible with regard to the horizontal resolution
and in 2008 the EMEP domain was extended to include EECCA countries.
The current version, and the provided gridded input data,
have a horizontal resolution of 50x50 km2 (at 60°N) and are defined
on a polar stereographic projection.

**IMPORTANT:**
The input data available in this ftp site should be appropriately acknowledged
when used for model runs.
If nothing else is specified according to references further in this chapter,
please acknowledge EMEP/MSC-W in any use of these data.

[CLRTAP]:   http://www.unece.org/env/lrtap/welcome.html
[publ2016]: http://emep.int/publ/emep2016_publications.html
[GPLv3]:    http://www.gnu.org/copyleft/gpl.html
[netCDF_CF]:http://www.unidata.ucar.edu/software/netcdf/conventions.html
[guide]:    http://github.com/metno/emep-ctm/releases/download/rv4_8/User_Guide.pdf
[rel410]:   http://github.com/metno/emep-ctm/releases/tag/rv4_10
[rel48]:    http://github.com/metno/emep-ctm/releases/tag/rv4_8
[rel45]:    http://github.com/metno/emep-ctm/releases/tag/rv4_5
[rel44]:    http://github.com/metno/emep-ctm/releases/tag/rv4_4
[rel43]:    http://github.com/metno/emep-ctm/releases/tag/rv4_3
[readme43]: http://github.com/metno/emep-ctm/releases/download/rv4_3/README_rv4_3special.txt
[rel40]:    http://github.com/metno/emep-ctm/releases/tag/rv4_0
[rel201106]:http://github.com/metno/emep-ctm/releases/tag/v201106
[rel30]:    http://github.com/metno/emep-ctm/releases/tag/rv3

**Previous releases (YYYYMM - date of release):**
* [OpenSource rv4.8 (201510)][rel48].
* [OpenSource rv4.5 (201409)][rel45].
* [OpenSource rv4.4 (201309)][rel44].  
* [OpenSource rv4.3 (201304)][rel43] *read the file [README_rv4_3special.txt][readme43]*.
* [OpenSource rv4.0 (201209)][rel40].
* [OpenSource v.2011-06 (201108)][rel201106].
* [OpenSource rv3 (200802)][rel30].

Please read the updated EMEP MSC-W model User Guide available
for each releaseon respective ftp site.

The EMEP/MSC-W model is designed to calculate air concentrations
and deposition fields for major acidifying and eutrophying pollutants,
photo-oxidants and particulate matter.
The present release of the model is intended to:
* facilitate insight in the model assumptions, the parametrisation used,
  the requirements for input data and the actual model code
* encourage dialogue and collaboration with the modelling community
* allow individual model runs and insight on how to run different scenarios

Questions on the EMEP model can be addressed to <emep.mscw@met.no>.
Support to the user community will be developed here with your contribution.
Please let us know what your needs for information are.

## Documentation
The EMEP/MSC-W model is a chemical transport model developed at the
Meteorological Synthesizing Centre - West (MSC-W)
at the Norwegian Meteorological Institute (met.no).
The EMEP model is a limited-area, terrain following sigma coordinate model
designed to calculate air concentration and deposition fields for

* acidifying and eutrophying compounds (S, N)
* ground level ozone (O3)
* particulate matter (PM2.5, PM10).

as well as their long-range transport and fluxes across national boundaries
(Transboundary air pollution).
A history of the development of the EMEP model can be found at
[http://www.emep.int/models][mscwmodels].

[mscwmodels]: http://www.emep.int/mscw/models.html#mscwmodels

### Model Description

The EMEP MSC-W chemical transport model -- technical description
D. Simpson, A. Benedictow, H. Berge, R. Bergstrõm, L. D. Emberson, H. Fagerli,
C. R. Flechard, G. D. Hayman, M. Gauss, J. E. Jonson, M. E. Jenkin, A. Nyíri,
C. Richter, V. S. Semeena, S. Tsyro, J.-P. Tuovinen, Á. Valdebenito, and P. Wind
Atmos. Chem. Phys., 12, 7825-7865, 2012
http://www.atmos-chem-phys.net/12/7825/2012/acp-12-7825-2012.html

### Computer requirements
To compile the EMEP model you need:

* Fortran 95 compiler
* NetCDF Library (>4.1.3)
* MPI Library (>1.0)

It is necessary to compile with double precision reals (8 bytes reals).
The program has been used on computers ranging from a Linux laptop
to supercomputers (Itanium2 cluster, Intel Xeon cluster, Cray XT4, IBM power5+).
It is compatible with all compilers tested so far: Intel, PGI, gfortran, XL fortran.
A Makefile is included, the path to netcdf (INCL and LLIB) has to be adapted
to your machine, and to the fortran compiler (F90) and flags (F90FLAGS)
to the compiler you are using.

The code has been tested with 1 to 1024 CPUs, and scales well (for large grids).
If only one CPU is used 1-2 GB memory is required.
If more than one, for example 64 CPUs are used, 200 MB of memory per CPU is enough
(in the case of a 132 X 159 grid size).
For runs on more than 32 CPUs, a fast interconnect is recommended
(infiniband for example), for smaller runs, gigabit ethernet is sufficient.
It takes ~3.5 hrs on 64*Xeon X5355 (2.66GHz) for a 1-year simulation.

When downloading input data in order to do a "base run" please make sure that there
are 35 Gb disc space available, especially due to large meteorology input files.
The model can be run for shorter periods, users can download meteorology for
only the period they are interested in, plus one day.

## Verification
The EMEP/MSC-W model is validated and reported to the
Cooperative Programme for Monitoring and Evaluation of the
Long-range Transmission for Air Pollutants in Europe ([EMEP][]) each year
by the EMEP/MSC-W group.
The reports can be found under the following links:
* [Status Reports][]
* [Country Reports][]

[EMEP]:           http://www.emep.int/
[Status Reports]: http://www.emep.int/publ/common_publications.html
[Country Reports]:http://www.emep.int/mscw/mscw_datanotes.html
