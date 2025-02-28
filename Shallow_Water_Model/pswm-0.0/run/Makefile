# Makefile for running PSWM:  Periodic Shallow-Water Model

# To run the model:  
#
# Create a file <runid>.dat containing the parameters for the run,
# then type "make RUNID=<runid>" where <runid> is your runid.
# If RUNID is not specified on the command line, it is taken from the file
# runid.txt; if that file doesn't exist, runid=cstest is used.
#
# Note that running the model sets the runid in runid.txt to that runid.
#
# See also README.txt for further information

BINDIR = ../bin
MODEL = pswm

# set the runid
ifndef ($(RUNID))
   RUNID := $(shell cat runid.txt)
   ifeq ($(strip $(RUNID)),)
      RUNID := pstest
   endif
endif


# To run the model for the run specified by RUNID
.PHONY: run
run : 
	-rm runid.txt
	echo $(RUNID) > runid.txt
	time $(BINDIR)/$(MODEL) $(RUNID) > output &

# To test the model
.PHONY: test
test : 
	-rm runid.txt pstest.out pstest_*.t[0-9]* output
	echo pstest > runid.txt
	@echo Running test run pstest--may take a minute...
	time $(BINDIR)/$(MODEL) pstest > output
	@echo ...done. Differences with archived test results:
	diff pstest.out pstest-output/pstest.out | more

# To clean up one run:  erase output from the run specified by RUNID
.PHONY: clean
clean:
	-rm $(RUNID).out $(RUNID)_*.t[0-9]* $(RUNID)_*.[1-9]* output

# To clean up all runs:  erase all output from all runs in this directory
.PHONY: clobber
clobber:
	-rm *.out *.t[0-9]* *.[1-9]* output

.PHONY : help
help: 
	@echo 'Running PSWM:  Periodic Shallow Water Model'
	@echo 
	@echo 'make                   to run the model with runid from runid.txt'
	@echo 'make RUNID=runid       to run the model with a specified runid'
	@echo 'make clean             to clean up from the last run'
	@echo 'make clean RUNID-runid to clean up from a specified run'
	@echo 'make clobber           to clean up from all runs'
	@echo 'make help              to see this list again'
	@echo 
	@echo '(see also the file README.txt)'

