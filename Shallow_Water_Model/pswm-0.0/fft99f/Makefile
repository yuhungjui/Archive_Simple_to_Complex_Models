# Standard Makefile for Library Packages (Fortran)
# Requires GNU make and GNU ar

#------------------------------------------------------------------------
# Usage:
#    list library source files in SRC
#    list additional files in ETC
#    specify destination directory for library in libdir
#    include "PRECISION=single" on make command line for single precision
#    use "make help" to list principal targets and actions
# Make variables used:
#    FFLAGS	flags passed to fortran compiler
#    LDFLAGS	flags passed to loader
#------------------------------------------------------------------------

# Author:
#
# Scott R. Fulton
# Department of Mathematics and Computer Science
# Clarkson University, Potsdam, NY    13699-5815
# fulton@clarkson.edu   www.clarkson.edu/~fulton

# Revision history:
# 15 Sep 2005:  Started this from a previous version 
# 21 Sep 2005:  Added support for multiple precision
# 25 Jun 2007:  Revised this copy for use with pswm
# 03 Aug 2007:  Included Makefile-fortran and Makefile-host

#------------------------------------------------------------------------
# preliminaries
#------------------------------------------------------------------------

# specify the shell (for safety)
SHELL	= /bin/sh

# specify some programs used
TAR	= tar
SED	= sed

# specify the date string used in dist files
#DATE	= `date '+%m%d%y'`
DATE	= `date '+%d%b%Y'`

# specify the precision [can override with PRECISION=single on command line] 
PRECISION	= double

#------------------------------------------------------------------------
# Details specific to this library
#------------------------------------------------------------------------

# set name of the library to that of the directory
CWD	= $(shell pwd)
package	= $(basename $(notdir $(CWD)))
LIB	= lib$(package).a

# specify where the completed library should go (could replace with ?)
libdir	= ../lib

# source files for this library
SRC	= fax.F fft991.F fft99a.F fft99b.F fft99.F fftfax.F fftrig.F vpassm.F

# corresponding object files (library members)
SRC_f	= $(SRC:.F=.f)
OBJ	= $(SRC:.F=.o)

# binary (executable) files (test programs)
BIN	= 

# extra files to be included in the distribution (package.tar.gz) file
ETC	= README.txt Log.txt Makefile
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# System-dependent information
#------------------------------------------------------------------------

include ../bin/Makefile-host

#------------------------------------------------------------------------
# Fortran rules
#------------------------------------------------------------------------

include ../bin/Makefile-fortran

#------------------------------------------------------------------------
# replace GNU make implicit rules to use sed instead of cpp for .F files
#------------------------------------------------------------------------

# sed script to convert .F source to single or double precision
ifeq ($(PRECISION),single)
   SED_SCRIPT = -e 's/EorD/e/g' -e 's/FLOAT/real/g' \
		-e 's/SINGLE_OR_DOUBLE/real            /g'
else
   SED_SCRIPT = -e 's/EorD/d/g' -e 's/FLOAT/dble/g' \
		-e 's/SINGLE_OR_DOUBLE/double precision/g'
endif

#  compile from .F source [preprocess using sed to set the precision]
%.o: %.F
	$(SED) $(SED_SCRIPT) $*.F > $*.f
	$(COMPILE.f) $(OUTPUT_OPTION) $*.f
	$(RM) $*.f

#  preprocess .F source to get .f source
%.f: %.F
	$(SED) $(SED_SCRIPT) $*.F > $*.f

#  compile and link from .F source [preprocess using sed to set the precision]
%: %.F
	$(SED) $(SED_SCRIPT) $*.F > $*.f
	$(LINK.f) $(OUTPUT_OPTION) $*.f
	$(RM) $*.f

#------------------------------------------------------------------------
# default target:  build the library
#------------------------------------------------------------------------

.PHONY : all
all : $(LIB)

#------------------------------------------------------------------------
# Additional targets specific to this directory:
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# build the library (uses implicit rules in GNU make)
#------------------------------------------------------------------------

$(LIB) : $(LIB)($(OBJ)) 
	@echo made $@

#------------------------------------------------------------------------
# install the library
#------------------------------------------------------------------------

.PHONY : install
install : $(LIB)
	mv $(LIB) $(libdir)
	@echo installed $@ in $(libdir)

#------------------------------------------------------------------------
# uninstall the library
#------------------------------------------------------------------------

.PHONY : uninstall
uninstall : clean
	-rm $(libdir)/$(LIB)

#------------------------------------------------------------------------
# clean up from building the library
#------------------------------------------------------------------------

.PHONY : clean
clean :
	-rm *.o $(LIB)

.PHONY : clobber
clobber : clean
	-rm $(package).tar.gz $(BIN)

#------------------------------------------------------------------------
# make a distribution tar file
#------------------------------------------------------------------------

.PHONY : dist
dist :
	-$(TAR) -czf $(package).tar.gz $(SRC) $(ETC)
	@echo made $(package).tar.gz

$(package).f : $(SRC_f)
	cat $(SRC_f) > $@

#------------------------------------------------------------------------
# Test programs
#------------------------------------------------------------------------

.PHONY : test
test : $(BIN)
	@echo running test programs:  $(BIN)
	@for name in $(BIN); do echo; \
	echo running test program $$name; ./$$name; done

#------------------------------------------------------------------------
# debugging
#------------------------------------------------------------------------

.PHONY : debug
debug :
	@echo libdir=$(libdir)
	@echo package=$(package)
	@echo date=$(DATE)
	@echo SRC=$(SRC)
	@echo ETC=$(ETC)
	@echo OBJ=$(OBJ)
	@echo LIB=$(LIB)
	@echo BIN=$(BIN)
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# miscellaneous housekeeping chores
#------------------------------------------------------------------------

# check the version of make being used
version :
	@echo Make version $(MAKE_VERSION)
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# display information on using this Makefile
#------------------------------------------------------------------------

.PHONY : help
help :
	@echo Library Package $(package)
	@echo 
	@echo 'make                 to compile the package'
	@echo 'make test            to run all test programs'
	@echo 'make install         to compile and install'
	@echo 'make clean           to clean up files from compilation'
	@echo 'make dist            to create a distribution tar file'
	@echo 'make clobber         to clean up all files created'
	@echo 'make uninstall       to uninstall the package'
	@echo 'make help            to get this help'

