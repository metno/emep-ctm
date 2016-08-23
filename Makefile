#
#
PROG =	Unimod
###################################################

include Makefile.SRCS

###################################################
LIBS = -lnetcdf -lnetcdff
INCL = -I/global/apps/netcdf/4.2.1.1/intel/13.0/include
LLIB = -L/global/apps/netcdf/4.2.1.1/intel/13.0/lib

#LIBS = -lnetcdf
#INCL = -I/global/apps/netcdf/3.6.2/include
#LLIB = -L/global/apps/netcdf/3.6.2/lib

F90 = mpif90

#GNU gfortran compiler (version 4.4.3 or later) 
F90FLAGS =  -ffree-line-length-none -fdefault-real-8   -O3

#Intel ifort compiler
F90FLAGS =  -shared-intel -r8  -recursive -O3

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

