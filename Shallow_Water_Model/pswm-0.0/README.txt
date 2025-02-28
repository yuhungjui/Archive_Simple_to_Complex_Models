File README.txt for directory /home/fulton/Models/pswm

Periodic Shallow Water Model

Documentation:
-------------
For an overview of the model and code organization, see ./doc/pswm-doc.txt.
For the details of equation forms and transforms, see ./doc/pswm-doc.pdf.
For model details see comments in the source code (start with the modules
./src/pswm_pars.f90 and ./src/pswm_data.f90).

Directories:
-----------
src	source code for the model
lib	compiled library packages
bin	directory for compiling the model
run	directory for running the model
plt	Matlab functions for reading and plotting the model output
doc	model documentation
fft99f	source code for Temperton's FFT package
test	directory for testing parts of the code
bugs	directory for archiving unresolved bugs

Files:
-----
Log.txt 	log of work on the model
README.txt	this file

Notes:  
-----
Not all directories and files listed above are in the distribution.
This directory is under version control by Mercurial (hg).

Scott R. Fulton
Department of Mathematics
Clarkson University, Potsdam, NY    13699-5815
fulton@clarkson.edu   www.clarkson.edu/~fulton
phone:   315-268-2379       FAX:  315-268-2371
----------------------------------------------
