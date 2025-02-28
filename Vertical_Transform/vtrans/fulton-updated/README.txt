File README.txt for directory vtrans:  Vertical Normal Modes Transforms
Last updated:
Fri May 29 16:31:39  EDT  2009

This directory contains software for vertical transforms as described in 

     Fulton, S. R. and W. H. Schubert, 1985:  Vertical Normal Mode Transforms:
     Theory and Application.  Monthly Weather Review, vol. 113, pp. 647-658.

The source files are:

   vtmain.f     main program
   vtrans.f     the vertical transform package
   cgamma.f     vertical structure functions for constant Gamma
   rsg.f        eispack routine for generalized eigenvalue problems
   second.f     routine for tracking execution times
   gaussl.f     library code for Gauss-Legendre quadrature
   testgl.f     simple test program for gaussl.f

All of these have been minimally updated so they should probably work on UNIX
systems.  Use of the various routines is fully explained in their comments.

One important caveat:  this is old Fortran77 code, written before Fortran90 
and before I knew anything about UNIX (or programming or numerical methods 
or ...) so don't expect too much.  Note in particular that the main program 
vtmain.f is quite a mess.

One more:  the code vtmain calls a function "adquad", which was an adaptive
quadrature (i.e., numerical integration) routine available from NCAR in the
1980s.  I don't have a copy any more, so this should be replaced by something
similar (perhaps from www.netlib.org).

Other files in this directory include:

   vttest.dat   input (data) file used to test the code
   mwr85.dat    input (data) file used for 1985 paper
   Makefile     makefile for (GNU make) for compiling, packaging, etc.
   Notes.txt    some brief comments on changes made
   README.txt   this file

The original code (in the subdirectory code1985) consists of:

   vtrans.for     Vertical transform package
   vtrans.dri     Driver program for use with data cards/files
   vtrans.tes     Sample job file (used to produce a test run)
   vtrans.job     Sample job file (used to produce results for paper above)
   cgamma.for     Vertical structure functions for constant Gamma
   cgamma.fod     Double precision version of  cgamma.for

Other subdirectories (possibly not present in tar or zip files) are:

   Dostalek/      files related to work by Jack Doestlek
   Archive/       copies of older versions
   Misc/          some other library routines not currently used

The code should run (I hope); I will try to help out if there are problems, 
but do not have time to fully support this package. 

Scott R. Fulton
Clarkson University
http://www.clarkson.edu/~fulton
fulton@clarkson.edu
