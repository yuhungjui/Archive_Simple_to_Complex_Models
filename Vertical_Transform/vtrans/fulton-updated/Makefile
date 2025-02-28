# Makefile for vtrans (vertical transforms)
# Written by Scott R. Fulton

NAME	= vtrans
SRC	= vtmain.f vtrans.f cgamma.f gaussl.f rsg.f second.f adquad.f
EXTRA	= Makefile *.f *.txt *.dat
OBJ	= $(SRC:.f=.o)
LIBS	=

FC	= gfortran
FFLAGS	= -O2
LDFLAGS	= 

# run make quietly...
.SILENT :

# create executable file:  simply say "make"
$(NAME) : $(OBJ)
	$(LINK.f) $(OBJ) $(LIBS) -o $(NAME)

# test program for Gauss-Legendre quadrature routine  gaussl
testgl : testgl.o gaussl.o
	$(LINK.f) testgl.o gaussl.o -o testgl

# remove object files
clean : 
	$(RM) *.o

# remove executable and object files
Clean : 
	$(RM) *.o $(NAME)
