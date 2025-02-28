# Makefile for library package:  vtrans (vertical transforms)
# Written by Scott R. Fulton
# Probably assumes GNU make (since that's what I use)

NAME	= vtrans

SRC	= vtmain.f vtrans.f cgamma.f gaussl.f rsg.f second.f

EXTRA	= Makefile *.f *.txt *.dat code1985/*

OBJ	= $(SRC:.f=.o)

LDFLAGS		= 

#LIBS	= -L$(HOME)/lib -lsrf

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

#----------------------------- zip files -------------------------------

ZIPFILES =      $(SRC) $(EXTRA)

NOZIP    =      \*.aux \*.log \*.dvi \*.b \*.bak \*.swp \*.zip \
		Archive Archive/* \*.tar \*.gz \*.tgz

# make a zip file of the complete source code:  make zip
# to include the date in the zip file name:     make zipd

ZIP	= zip

.PHONY : zip zipd

zip : $(ZIPFILES)
	$(ZIP) -u $(NAME).zip $(ZIPFILES) -x $(NOZIP)
	@echo made $(NAME).zip

zipd : $(ZIPFILES)
	$(ZIP) -u $(NAME)`date '+%m%d%y'`.zip $(ZIPFILES) -x $(NOZIP)
	@echo made $(NAME)`date '+%m%d%y'`.zip

#----------------------------- tar files -------------------------------

TARFILES =      $(SRC) $(EXTRA)

NOZIP    =      \*.aux \*.log \*.dvi \*.b \*.bak \*.swp \*.zip \
		Archive Archive/* \*.tar \*.gz \*.tgz

# make a gzipped tar file of the complete source code:  make tar
# to include the date in the gzipped tar file name:     make tard

TAR	= tar

.PHONY : tar tard

tar : $(TARFILES)
	$(TAR) -cvzf $(NAME).tar.gz $(TARFILES)
	@echo made $(NAME).tar.gz

tard : tar
	mv $(NAME).tar.gz $(NAME)`date '+%m%d%y'`.tar.gz
	@echo made $(NAME)`date '+%m%d%y'`.tar.gz

#-----------------------------------------------------------------------

# make list of complete source code (for printing):  make NAME.slist

slist $(NAME).slist : $(SRC) $(EXTRA)
	fjoin $(SRC) $(EXTRA) > $(NAME).slist

#---------------------------------------------------------------------
# macros needed for SUN workstations running under AFS
#LINK.g		= /afs/cu/software/Apps/ncar4.1/@sys/bin/ncargf77
LINK.g		= ncargf77
GLIB		= -L$(HOME)/lib -lplotf 
#FFLAGS	= -u -fast
#FFLAGS	= -u -fast -xO5
#FFLAGS	= -u -fast -xO5 -xchip=ultra 
FFLAGS	= -u -fast -xO5 -xarch=v8plusa
#FFLAGS	= -u -fast -xO5 -xarch=v8plusa -xchip=ultra 
#---------------------------------------------------------------------
# macros needed for linux PC
#LINK.g		= ncargf77
#GLIB		= -L$(HOME)/lib -lplotf 
#---------------------------------------------------------------------
# macros needed for HP workstations
#COMPILE.f	= f77 $(FFLAGS) -c
#LINK.f		= f77 $(LDFLAGS)
#RM		= rm
#LINK.g		= ncargf77
#GLIB		= -lplotf 
#---------------------------------------------------------------------

