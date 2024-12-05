#
#
PROG =	emepctm
###################################################

#include Makefile.SRCS

###################################################

# prefered netCDF 4.2.1.1 or later
LIBS = -lnetcdf -lnetcdff
#explicit pathes needed only if nf-config does not work
INCL = -I/my/path/to/include
LLIB = -L/my/path/to/lib

# options using nf-config utility (older versions used nc-config)
INCL = -I$(shell nf-config --includedir)
LLIB = -L$(shell nf-config --flibs)

F90 = mpif90

# GNU gfortran compiler (tested for version 8.5.0)
F90FLAGS = -fdefault-real-8  -ffixed-line-length-none -ffree-line-length-none -Wno-error=line-truncation -O3 -g
# GNU gfortran compiler (tested for version 12.2.0 and 11.3.0)
F90FLAGS = -fdefault-real-8 -fallow-argument-mismatch  -ffixed-line-length-none -ffree-line-length-none -Wno-error=line-truncation -O3 -g
#DEBUG flag
#F90FLAGS += -Wall -fbacktrace -fbounds-check -fimplicit-none -pedantic

# Intel ifort compiler (comment out if gfortran used)
F90FLAGS = -g -r8 -IPF_fp_relaxed -assume noold_maxminloc -O2 -march=core-avx2

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
