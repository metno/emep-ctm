
### NEW in this RELEASE ###

This tar file of Opensource code rv4.3 contains the model code, the runscript, 
namelist file called 'config_emep.nml' and the GPL licence. 'config_emep.nml' 
is a new feature of the model and the necessary flags and switches are defined 
here, whereas it was defined in the module 'ModelConstants_ml.f90' in the 
previous versions.


### Important MESSAGES about model code rv4.3 ###

An error was found recently with the way soil-NOx emissions are working. 
These emissions depend upon a pre-calculated field of annual N-depositions. 
To save having to re-calculate that field every time a new future scenario
was introduced, i.e., a scaling system which compares the future EU emissions 
with current-day ones, and just scaled accordingly. The problem with this is 
that if changes are made in emissions from just one EU country, the whole 
soil-NOx emission field changes, so we can get emission changes in e.g. the 
Po valley. So care should be taken while using this code for major production 
runs such as Source Receptor calculations.

Also please note that the code does not work well with 'gfortran'. All these 
will be fixed soon and a new Open Source code will be released after EMEP 
Reporting this year, which is in June 2013. In the mean time, a revised version 
of the code will be available on request for 'gfortran' and Source Receptor
purposes.
