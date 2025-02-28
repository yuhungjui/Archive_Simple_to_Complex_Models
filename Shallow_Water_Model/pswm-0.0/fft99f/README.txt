File README.txt for directory /home/fulton/Models/pswm/fft99f

FFT99F:  All-Fortran version of Temperton's FFT routines (see also ecmfft)

Note:  This is my repackaging of source code available from the web.

Documentation:
-------------
Primary documentation is comments in file fft99.F.

Contents:
--------
*.F		source code (variable precision)
Makefile	makefile for compiling and installing library
Log.txt		log of changes, etc.
README.txt	this file

Distribution:  gzipped tar file fft99f.tar.gz ("make dist" to generate)

To install:
	tar -xzf fft99f.tar.gz
	make install

Note:
----
This is a copy of the primary source code for fft99f in ~/lib/fft99f.
All extraneous files (previous versions, output, etc.) have been removed
and the Makefile has been adjusted to install fft99f where needed for pswm.

Scott R. Fulton
Department of Mathematics and Computer Science
Clarkson University, Potsdam, NY    13699-5815
fulton@clarkson.edu   www.clarkson.edu/~fulton
phone:   315-268-2379       FAX:  315-268-2371
----------------------------------------------

