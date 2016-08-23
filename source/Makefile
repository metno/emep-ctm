#
#
PROG =	Unimod
###################################################

include Makefile.SRCS

###################################################

LIBS = -lnetcdf -lnetcdff
INCL = -I/global/apps/netcdf/4.2.1.1/intel/13.0/include
LLIB = -L/global/apps/netcdf/4.2.1.1/intel/13.0/lib


F90 = mpif90

#Intel ifort compiler
F90FLAGS =  -shared-intel -r8  -recursive -O3


###GNU gfortran compiler (version 4.4.3 or later) 
#LIBS = -lnetcdf -lnetcdff
#INCL = -I$(NETCDF_ROOT)/include -I$(filter /global/apps/openmpi%,$(subst :, ,$(CPATH)))
#LLIB = -L$(NETCDF_ROOT)/lib

#F90FLAGS =  -ffree-line-length-none -fdefault-real-8   -O3


###################################################


LDFLAGS = $(F90FLAGS) $(LLIB) -o $(PROG) $(FOBJ) $(INCL) $(LIBS)


.SUFFIXES: $(SUFFIXES)  .f90

.f90.o:
	$(F90) $(F90FLAGS) $(INCL) -c $<


all:  $(PROG)


# Include the dependency-list (created by makedepf90)
include dependencies

$(PROG): $(FOBJ)
	 $(F90) $(LDFLAGS)
#

clean: diskclean

diskclean:
	rm -f $(PROG) *.o *.mod

##########################################################

