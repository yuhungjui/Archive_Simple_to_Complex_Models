SHELL=/bin/tcsh
SRCS=	gsw.f baleqdt.f bogusuv.f cons.f getdta.f hdiffu.f \
	lptime.f lrhs.f nrhs.f ntulib1.f outflds.f \
	prepare.f rdy2go.f rk4time.f track.f fftx.f getvor1.f \
        qmaxn3.f uven0.f uven80.f getvor2.f

OBJS=	gsw.o baleqdt.o bogusuv.o cons.o getdta.o hdiffu.o \
	lptime.o lrhs.o nrhs.o ntulib1.o outflds.o \
	prepare.o rdy2go.o rk4time.o track.o fftx.o getvor1.o \
        qmaxn3.o uven0.o uven80.o getvor2.o

FC =		pgf77 
CMD =		gsw.exe
LIBS =

FFLAGS	 = -r8 -O
                                                                  
LDFLAGS	= -r8 -O

all:		$(CMD)

$(CMD):		$(OBJS)
		$(FC) -o $(@) $(LDFLAGS) $(OBJS) $(LIBS)
                                                                
clean:
	-rm -f $(OBJS)

clobber:	clean
	-rm -f $(CMD)
                
void:	clobber
	-rm -f $(SRCS)
