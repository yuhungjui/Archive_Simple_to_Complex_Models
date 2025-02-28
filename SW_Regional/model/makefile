SHELL=/bin/sh
#
SRCS=	block.f boundary.f bogusuv.f getum.f nbal.f \
	nbalm.f outflds.f qmaxn3.f rcons.f regrhs.f gettopo.f \
	rgetrdy.f rhdiffu.f rintgrt.f rsw.f reglib.f fftx.f \
	bogusbanduv.f bogusuv2.f bogusty.f getvorm1.f
#
OBJS=	block.o boundary.o bogusuv.o getum.o nbal.o \
	nbalm.o outflds.o qmaxn3.o rcons.o regrhs.o gettopo.o \
	rgetrdy.o rhdiffu.o rintgrt.o rsw.o reglib.o fftx.o \
	bogusbanduv.o bogusuv2.o bogusty.o getvorm1.o

# Tunable parameters
#
# FC		Name of the fortran compiling system to use
# LDFLAGS	Flags to the loader
# LIBS		List of libraries
# CMD		Name of the executable
# PROFLIB	Library needed for profiling
#
#FC =		cf77
FC =		pgf77
LDFLAGS =	-O2
LIBS =		
CMD =		regtopo.exe
PROFLIB =	-lprof

# To perform the default compilation, use the first line
# To compile with flowtracing turned on, use the second line
# To compile giving profile additonal information, use the third line
# WARNING:  SIMULTANEOUSLY PROFILING AND FLOWTRACING IS NOT RECOMMENDED 
#FFLAGS = -I$(INCLUDE) -Ad -Sw
#FFLAGS = -I$(INCLUDE) -Ad -Sw -Of -Wv,-Of
##FFLAGS = -I$(INCLUDE) -Ad -Sw -Wv,-Of
#FFLAGS =	 -Zv
#FFLAGS =	 -F
#FFLAGS =	 -Wf"-ez"

# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build a.out.
#
all:		$(CMD)

$(CMD):		$(OBJS)
	$(FC) $(LDFLAGS) -o $(@) $(OBJS) $(LIBS)

# Make the profiled version of the command and call it a.out.prof
#
$(CMD).prof:	$(OBJS)
	$(FC) $(LDFLAGS) -o $(@) $(OBJS) $(PROFLIB) $(LIBS)

clean:
	-rm -f $(OBJS)

clobber:	clean
	-rm -f $(CMD) $(CMD).prof

void:	clobber
	-rm -f $(SRCS) makefile
