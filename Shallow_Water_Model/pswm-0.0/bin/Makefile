# Makefile for /home/fulton/Models/pswm/bin
#
# Requires GNU version of make (probably at least version 3.80)
# Type "make help" for list of main targets and actions

#-----------------------------------------------------------------------------
# Author:
#
# Scott R. Fulton
# Department of Mathematics and Computer Science
# Clarkson University, Potsdam, NY    13699-5815
# fulton@clarkson.edu   www.clarkson.edu/~fulton

# Revision history:
# 30 Jul 2007:  Original version

#-----------------------------------------------------------------------------
# Utilities used
#-----------------------------------------------------------------------------

# specify the shell (for safety)
SHELL	= /bin/sh

#-----------------------------------------------------------------------------
# other Makefile variables
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Details specific to this model
#-----------------------------------------------------------------------------

MODEL	= pswm
BASEDIR	= ..

# specify the directories for source, libraries, and compiling
SRCDIR  = $(BASEDIR)/src
LIBDIR  = $(BASEDIR)/lib
BINDIR  = $(BASEDIR)/bin

# model name (used for executable code)
MODEL	= pswm
BIN	= $(MODEL)

# default target:  build the model
.PHONY : all
all : $(BIN)

# source files for this model
SRC =	pswm_main.f90 pswm_pars.f90 pswm_cons.f90 pswm_vars.f90 \
	pswm_setup.f90 pswm_terms.f90 pswm_output.f90 pswm_data.f90 \
	ss_ops.f90 sitpack_interface.f90 sitpack.f90 kinds.f90 

# corresponding object files
OBJ = $(SRC:.f90=.o)

# dependencies (generated using getdep -M)
kinds.o: kinds.f90
pswm_cons.o: kinds.o pswm_cons.f90
pswm_data.o: kinds.o pswm_pars.o pswm_cons.o ss_ops.o pswm_data.f90
pswm_main.o: kinds.o sitpack_interface.o sitpack.o pswm_pars.o pswm_vars.o pswm_setup.o pswm_data.o pswm_output.o pswm_main.f90
pswm_output.o: kinds.o pswm_pars.o pswm_cons.o pswm_vars.o pswm_terms.o ss_ops.o pswm_output.f90
pswm_pars.o: kinds.o pswm_pars.f90
pswm_setup.o: kinds.o pswm_pars.o pswm_cons.o pswm_vars.o pswm_terms.o sitpack_interface.o sitpack.o ss_ops.o pswm_setup.f90
pswm_terms.o: kinds.o pswm_pars.o pswm_cons.o ss_ops.o pswm_data.o pswm_terms.f90
pswm_vars.o: kinds.o pswm_vars.f90
read_icpars.o: read_icpars.f90
sitpack.o: kinds.o sitpack_interface.o sitpack.f90
sitpack_interface.o: kinds.o pswm_pars.o pswm_cons.o pswm_vars.o pswm_terms.o sitpack_interface.f90
ss_ops.o: kinds.o ss_ops.f90

# search path for other source files (here, just for sitpack.f90)
vpath %.f90 $(SRCDIR)

# libraries:  files, link names, and names
LIBFILES	= libfft99f.a
LIBLINKS	= $(LIBFILES:lib%.a=-l%)
LIBNAMES	= $(LIBFILES:lib%.a=%)

# library search path and library names for loader (used in LINK.* commands)
LOADLIBES	= -L$(LIBDIR)
LDLIBS		= $(LIBLINKS)

#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# System-dependent information (complier names, options, etc.)
#-----------------------------------------------------------------------------

include Makefile-host

#------------------------------------------------------------------------
# Rules and commands for processing Fortran
#------------------------------------------------------------------------

include Makefile-fortran

#-----------------------------------------------------------------------------
# build the model (uses implicit rules in GNU make)
#-----------------------------------------------------------------------------

$(MODEL) : $(OBJ) $(LIBFILES)
	$(LINK.f) $(OBJ) $(LOADLIBES) $(LDLIBS) -o $@
	@echo made $@


.PHONY : install
install : $(MODEL) clean

#-----------------------------------------------------------------------------
# build the libraries (GNU make static pattern rule)
#-----------------------------------------------------------------------------

$(LIBFILES) : lib%.a :
	cd $(BASEDIR)/$*; make install

#-----------------------------------------------------------------------------
# test the model
#-----------------------------------------------------------------------------

.PHONY : check
check : $(MODEL)
	@-rm output_new
	./$(MODEL) > output_new
	diff output output_new | more
	@echo 'Output of $(MODEL) (output_new) should match previous results (output)'.

#-----------------------------------------------------------------------------
# clean up from building the model
#-----------------------------------------------------------------------------

.PHONY : clean
clean :
	-rm *.o *.mod

.PHONY : uninstall clobber
uninstall clobber : clean
	-rm $(BIN)

#-----------------------------------------------------------------------------
# debugging (for this Makefile)
#-----------------------------------------------------------------------------

.PHONY : debug
debug :
	@echo MODEL=$(MODEL)
	@echo SRC=$(SRC)
	@echo OBJ=$(OBJ)
	@echo BIN=$(BIN)
	@echo SRCDIR=$(SRCDIR)
	@echo LIBDIR=$(LIBDIR)
	@echo BINDIR=$(BINDIR)
	@echo LIBFILES=$(LIBFILES)
	@echo LIBLINKS=$(LIBLINKS)
	@echo LIBNAMES=$(LIBNAMES)
	@echo LOADLIBES=$(LOADLIBES)
	@echo LDLIBS=$(LDLIBS)

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# miscellaneous housekeeping chores
#-----------------------------------------------------------------------------

# check the version of make being used
version :
	@echo Make version $(MAKE_VERSION)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# display information on using this Makefile
#-----------------------------------------------------------------------------

.PHONY : help
help :
	@echo Model $(MODEL)
	@echo 
	@echo 'make           compiles the model'
	@echo 'make check     runs the model and checks its output'
	@echo 'make clean     cleans up object files from compilation'
	@echo 'make clobber   cleans up object files and compiled model'
	@echo 'make version   shows version of "make" used (should be >=3.80)'
	@echo 'make help      produces this help'
