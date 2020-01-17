# Open Source EMEP MSC-W model
[![Documentation Status](https://readthedocs.org/projects/emep-ctm/badge/?version=latest)](http://emep-ctm.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/65725508.svg)](https://zenodo.org/badge/latestdoi/65725508)


The EMEP models have been instrumental to the development of
air quality policies in Europe since the late 1970s,
mainly through their support to the strategy work under the
[Convention on Long-range Transboundary Air Pollution][CLRTAP].
In the 1990s the EMEP models became also the reference tools for
atmospheric dispersion calculations as input to the Integrated Assessment Modelling,
which supports the development of air quality polices in the European Union.

[CLRTAP]:   http://www.unece.org/env/lrtap/welcome.html
[GPLv3]:    http://www.gnu.org/copyleft/gpl.html
[netCDF_CF]:http://www.unidata.ucar.edu/software/netcdf/conventions.html
[guide]:    http://emep-ctm.readthedocs.io/en/latest
[publ2019]: http://emep.int/publ/emep2019_publications.html
[rel433]:   http://github.com/metno/emep-ctm/releases/tag/rv4_33
[rel432]:   http://github.com/metno/emep-ctm/releases/tag/rv4_32
[rel417]:   http://github.com/metno/emep-ctm/releases/tag/rv4_17
[rel415]:   http://github.com/metno/emep-ctm/releases/tag/rv4_15
[rel410]:   http://github.com/metno/emep-ctm/releases/tag/rv4_10
[rel48]:    http://github.com/metno/emep-ctm/releases/tag/rv4_8
[rel45]:    http://github.com/metno/emep-ctm/releases/tag/rv4_5
[rel44]:    http://github.com/metno/emep-ctm/releases/tag/rv4_4
[rel43]:    http://github.com/metno/emep-ctm/releases/tag/rv4_3
[readme43]: http://github.com/metno/emep-ctm/releases/download/rv4_3/README_rv4_3special.txt
[rel40]:    http://github.com/metno/emep-ctm/releases/tag/rv4_0
[rel201106]:http://github.com/metno/emep-ctm/releases/tag/v201106
[rel30]:    http://github.com/metno/emep-ctm/releases/tag/rv3

The EMEP MSC-W model is designed to calculate air concentrations
and deposition fields for major acidifying and eutrophying pollutants,
photo-oxidants and particulate matter.
The Open Source releases are intended to:
* facilitate insight in the model assumptions, the parametrisation used,
  the requirements for input data and the actual model code
* encourage dialogue and collaboration with the modelling community
* allow individual model runs and insight on how to run different scenarios

## Releases

The latest Open Source EMEP MSC-W model version ([rv4.33][rel433])
will be used on the [EMEP status reporting of the year 2019][publ2019].
The source code, together with a set of input data,
an updated user guide and a full year model results for the year 2015,
under [GPL license v3][GPLv3].

You can download the source code and input data by following the instructions on the 
[OpenSource rv4.33 (201906)][rel433] release page.

Information on use of the model can be found in the [EMEP MSC-W model User Guide][guide]

#### Previous releases (YYYYMM - date of release)

* [OpenSource rv4.32 (201904)][rel432].
* [OpenSource rv4.17 (201802)][rel417].
* [OpenSource rv4.15 (201709)][rel415].
* [OpenSource rv4.10 (201609)][rel410].
* [OpenSource rv4.8 (201510)][rel48].
* [OpenSource rv4.5 (201409)][rel45].
* [OpenSource rv4.4 (201309)][rel44].  
* [OpenSource rv4.3 (201304)][rel43] *read the file [README_rv4_3special.txt][readme43]*.
* [OpenSource rv4.0 (201209)][rel40].
* [OpenSource v.2011-06 (201108)][rel201106].
* [OpenSource rv3 (200802)][rel30].

### Model domain and resolution

The EMEP MSC-W model is very flexible with regard to the horizontal resolution
and vertical resolutions. In 2008 the EMEP domain was extended to include EECCA countries.
In 2017 the vertical resolution increased from 20 to 34 model levels.
In the latest release, the provided gridded input and output data
are provided on 2 different model domains and resolutions:
- `EECCA` domain with a horizontal resolution of 50x50 km2 (at 60°N), 
  on polar stereographic projection, and 20 vertical levels;
- `EMEP01` domain with a 0.1x0.1 degrees on long-lat projection,
  and 34 vertical levels.

### Input data

The standard input files for the EMEP model are provided in two different formats,
netCDF or ASCII.
We presently follow as much as possible the [netCDF CF conventions][netCDF_CF]
for both input and output data.
More details about the input files are described in Chapter 2 of the
[EMEP MSC-W model User Guide][guide].

**IMPORTANT:**
The input data that accompanies each release model should be appropriately acknowledged
when used for model runs.
If nothing else is specified according to references further in this chapter,
please acknowledge EMEP/MSC-W in any use of these data.

### Model output

The main output files are in netCDF format
and some additional output files are in ASCII format.
Details about the output files are described in Chapter 3 of the User Guide.

## Documentation
The EMEP MSC-W model is a chemical transport model developed at the
Meteorological Synthesizing Centre - West (MSC-W)
at the Norwegian Meteorological Institute (met.no).
The EMEP model is a limited-area, terrain following hybrid coordinate model
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

## Running the model

How to submit a run is described in detail in Chapter 4
of the [EMEP MSC-W model User Guide][guide].

Questions on the EMEP model can be submitted to the [issues tracker].
Support to the user community will be developed here with your contribution.
Please let us know what your needs for information are.

[issues tracker]: https://github.com/metno/emep-ctm/issues

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

The EMEP MSC-W model is validated and reported to the
Cooperative Programme for Monitoring and Evaluation of the
Long-range Transmission for Air Pollutants in Europe ([EMEP][]) each year
by the EMEP/MSC-W group.
The reports can be found under the following links:
* [Status Reports][]
* [Country Reports][]

[EMEP]:           http://www.emep.int/
[Status Reports]: http://www.emep.int/publ/common_publications.html
[Country Reports]:http://www.emep.int/mscw/mscw_datanotes.html

## Citation

```bibtex
@article{emep-mscw_2012,
  author = {Simpson, D. and Benedictow, A. and Berge, H. and Bergstr\"om, R. and Emberson, L. D. and 
            Fagerli, H. and Flechard, C. R. and Hayman, G. D. and Gauss, M. and Jonson, J. E. and
            Jenkin, M. E. and Ny\'{\i}ri, A. and Richter, C. and Semeena, V. S. and Tsyro, S. and
            Tuovinen, J.-P. and Valdebenito, \'A. and Wind, P.},
  title = {The EMEP MSC-W chemical transport model &ndash; technical description},
  journal = {Atmospheric Chemistry and Physics},
  volume = {12},
  year = {2012},
  number = {16},
  pages = {7825--7865},
  doi = {10.5194/acp-12-7825-2012},
  url = {https://www.atmos-chem-phys.net/12/7825/2012/}
}

@misc{emep-ctm_rv4.33,
  author       = {EMEP MSC-W},
  title        = {Open Source EMEP/MSC-W model rv4.33 (201906)},
  month        = jun,
  year         = 2019,
  doi          = {10.5281/zenodo.3265912},
  url          = {https://doi.org/10.5281/zenodo.3265912}
}

@misc{emep-ctm_rv4.17,
  author       = {EMEP MSC-W},
  title        = {Open Source EMEP/MSC-W model rv4.17 (201802)},
  month        = feb,
  year         = 2018,
  doi          = {10.5281/zenodo.3355023},
  url          = {https://doi.org/10.5281/zenodo.3355023}
}

@misc{emep-ctm_rv4.15,
  author       = {EMEP MSC-W},
  title        = {Open Source EMEP/MSC-W model rv4.15 (201709)},
  month        = sep,
  year         = 2017,
  doi          = {10.5281/zenodo.3355041},
  url          = {https://doi.org/10.5281/zenodo.3355041}
}

@misc{emep-ctm_rv4.10,
  author       = {EMEP MSC-W},
  title        = {Open Source EMEP/MSC-W model rv4.10 (201609)},
  month        = sep,
  year         = 2016,
  doi          = {10.5281/zenodo.3355083},
  url          = {https://doi.org/10.5281/zenodo.3355083}
}
```
