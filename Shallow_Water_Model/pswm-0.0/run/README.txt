File README.txt for directory /home/fulton/Models/pswm/run

Directory for running PSWM:  Periodic Shallow Water Model

Overview:
--------
To run the model (using the supplied main program ../src/pswm_main.f90):

1. Create an input file with name <runid>.dat (where <runid> is the runid).
For example, there is a file here called pstest.dat for the runid "pstest".

2. Run the model by one of the following:
(a) Type "make RUNID=runid" (for example:  "make RUNID=pstest")
(b) Create a file "runid.txt" containing your runid type "make"
(c) Type "../bin/pswm runid" (for example:  "../bin/pswm pstest") 

3. To clean up after running the model do one of the following:
(a) Type "make clean RUNID=runid" to clean up from run runid
(b) Type "make clean" to clean up the latest run (runid as in runid.txt) 
(c) Type "make clobber" to clean up all model output (all runs)

Standard Runs:  There is a .dat file for each of the following:
-------------

  runid		description of run
  -----		------------------
  pstest	standard test run (results archived in ./pstest-output)
  (others to follow later)

Scott R. Fulton
Department of Mathematics
Clarkson University, Potsdam, NY    13699-5815
fulton@clarkson.edu   www.clarkson.edu/~fulton
phone:   315-268-2379       FAX:  315-268-2371
----------------------------------------------
