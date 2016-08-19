#
#
PROG =	Unimod
###################################################

include Makefile.SRCS

###################################################


F90_CONFORM_CHECK_ABORT=ON

###################################################

LIBS = -lmpi -lnetcdf
INCL = -I/home/u4/mifahik/netcdf/include
LLIB = -L/home/u4/mifahik/netcdf/lib64


F90 = f90

F90FLAGS = -64 -r8 -O2 -mieee-fp -ftz -OPT:IEEE_arithm=3:roundoff=3 -TARG:exc_min=0ZV $(INCL)
#F90FLAGS = -64 -r8 -O3 -OPT:IEEE_arithm=3:roundoff=3 -TARG:exc_min=0ZV $(INCL)
#F90FLAGS = -64 -r8 -g -C -DEBUG:trap_uninitialized=ON:verbose_runtime=ON -TARG:exc_min=0ZV $(INCL)
#F90FLAGS = -64 -r8 -g -C -DEBUG:trap_uninitialized=ON:verbose_runtime=ON -DEBUG:conform_check=ON -DEBUG:subscript_check:verbose_runtime=ON -DEBUG:fullwarn=ON -TARG:exc_min=0ZV $(INCL)

LDFLAGS = -64 -r8 -O2 -mieee-fp -ftz  -OPT:IEEE_arithm=3:roundoff=3 -TARG:exc_min=0ZV
#LDFLAGS = -64 -r8 -O 3 -OPT:IEEE_arithm=3:roundoff=3 -TARG:exc_min=0ZV
#LDFLAGS = -64 -r8 -C -g -DEBUG:trap_uninitialized=ON:verbose_runtime=ON -TARG:exc_min=0ZV
#LDFLAGS = -64 -r8 -g -C -DEBUG:trap_uninitialized=ON:verbose_runtime=ON -DEBUG:conform_check=ON -DEBUG:subscript_check:verbose_runtime=ON -DEBUG:fullwarn=ON -TARG:exc_min=0ZV
LD = f90


.SUFFIXES: $(SUFFIXES)  .f .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

.f.o:	
	$(F90) $(F90FLAGS) -c $<


# Include the dependency-list created by makedepf90 below
all:  $(PROG)

include .depend

#LLIB added, ds
depend .depend:
	/home/u4/mifahik/bin/makedepf90 $(SRCS) \
	-o $(PROG) \
	-l "$$(F90) $$(LDFLAGS) $$(LLIB) -o $$(PROG) $$(FOBJ) $$(INCL) $$(LIBS)" > .depend

clean:
	rm -f $(PROG) *.o *.mod .depend; \
	#touch .depend
#make depend


##########################################################

