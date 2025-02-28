LIBNAME = libverttrans_dble.a

OBJ     = vtrans.o gaussl.o rsg.o
          
LIBOBJ  = $(OBJ) 

FC     = pgf90
FFLAGS = -O2 

$(LIBNAME): $(LIBOBJ)
	ar r $@ $(LIBOBJ)
	ranlib $@

%.o: ../%.f
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -f *.o $(LIBNAME)
