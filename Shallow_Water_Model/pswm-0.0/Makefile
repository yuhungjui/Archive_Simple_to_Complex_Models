# Makefile for PSWM:  Periodic Shallow Water Model
# Requires GNU version of make, tar, and gzip

#------------------------------------------------------------------------
# Usage:
#    specify destination directory for library in libdir
#    use "make help" to list principal targets and actions
#------------------------------------------------------------------------

# Author:
#
# Scott R. Fulton
# Department of Mathematics and Computer Science
# Clarkson University, Potsdam, NY    13699-5815
# fulton@clarkson.edu   www.clarkson.edu/~fulton

# Revision history:
# 10/20/05:  Original version (for cswm)
# 12/07/05:  Added support for distribution and backup files
# 06/06/07:  Fixed a couple of trivial bugs
# 08/15/07:  Updated for use with pswm

#------------------------------------------------------------------------
# preliminaries
#------------------------------------------------------------------------

# specify the shell (for safety)
SHELL	= /bin/sh

# specify some programs used
TAR	= tar

# specify the date string used in dist and backup files
#DATE	= `date '+%m%d%y'`
#DATE	= `date '+%d%b%Y'`
DATE	= `date '+%F'`

#------------------------------------------------------------------------
# Directories
#------------------------------------------------------------------------
# set name of the model to that of the directory
MODEL	= pswm
CWD	= $(shell pwd)
BASEDIR	= $(CWD)
LIBDIR	= $(BASEDIR)/lib
BINDIR	= $(BASEDIR)/bin
RUNDIR	= $(BASEDIR)/run
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# default target:  print a summary of options and help
#------------------------------------------------------------------------

.PHONY : welcome
welcome help: 
	@echo 'PSWM:  Periodic Shallow Water Model'
	@echo 
	@echo 'For help see files README.txt and doc/pswm-doc.txt'
	@echo 
	@echo 'To use the model:'
	@echo 'make install    to compile the model'
	@echo 'make test       to test the model'
	@echo 'make uninstall  to clean up all files created'
	@echo 'make dist       to create a distribution tar file'
	@echo 'make backup     to create a backup tar file (all files)'
	@echo 'make help       to see this list again'

#------------------------------------------------------------------------
# install the model
#------------------------------------------------------------------------

.PHONY : install
install : 
	cd $(BINDIR); make install

#------------------------------------------------------------------------
# test the model
#------------------------------------------------------------------------

.PHONY : test
test : install
	cd $(RUNDIR); make test

#------------------------------------------------------------------------
# uninstall the model (and libraries)
#------------------------------------------------------------------------

.PHONY : uninstall-model
uninstall-model : 
	-cd $(RUNDIR); make clobber
	-cd $(BINDIR); make uninstall

.PHONY : uninstall
uninstall : uninstall-model
	-cd $(LIBDIR); rm lib*.a

#------------------------------------------------------------------------
# make a distribution file (with version number)
#------------------------------------------------------------------------

TAROPTIONS = --exclude='*.a' --exclude='*.b' --exclude='Log.txt' \
	     --exclude='*.gz' --exclude='*.mat' --exclude='.hg*' \
	     --exclude='test*' --exclude='bugs*'

.PHONY : dist
dist : uninstall-model
	-$(TAR) -czf $(MODEL)-$(DATE).tar.gz $(TAROPTIONS) *
	@echo made $(MODEL)-$(DATE).tar.gz

#------------------------------------------------------------------------
# make a backup file (with date)
#------------------------------------------------------------------------

BACKOPTIONS = --exclude='*.a' --exclude='pswm-*.gz' --exclude='*.mat'

.PHONY : backup
backup : uninstall-model
	-$(TAR) -czf $(MODEL)-$(DATE).tar.gz $(BACKOPTIONS) *
	@echo made $(MODEL)-$(DATE).tar.gz

#------------------------------------------------------------------------
# debugging
#------------------------------------------------------------------------

.PHONY : debug
debug :
	@echo MODEL=$(MODEL)
	@echo BASEDIR=$(BASEDIR)
	@echo LIBDIR=$(LIBDIR)
	@echo BINDIR=$(BINDIR)
	@echo RUNDIR=$(RUNDIR)
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# miscellaneous housekeeping chores
#------------------------------------------------------------------------

# check the version of make being used
version :
	@echo Make version $(MAKE_VERSION)
#------------------------------------------------------------------------
