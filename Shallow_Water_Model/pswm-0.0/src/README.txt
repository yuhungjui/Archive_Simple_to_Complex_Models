File README.txt for directory /home/fulton/Models/pswm/src

Source code for Periodic Shallow Water Model

Files
-----
pswm_main.f90		main program
pswm_pars.f90		model parameters (user-specifiable constants)
pswm_cons.f90		model constants  (not user-specifiable, constant)
pswm_vars.f90		model variables  (change during a run)
pswm_data.f90		data for a model run:  initial conditions and forcing
pswm_setup.f90		routines for initializing the model
pswm_terms.f90		primary code to compute the terms in the model
pswm_output.f90		output routines
ss_ops.f90		routines to compute operations in spectral space
sitpack_interface.f90	interface routines between the model and  sitpack
sitpack.f90		model-independent semi-implicit time integration package
kinds.f90		definitions of variable kinds used in the model
README.txt		this file

Scott R. Fulton
Department of Mathematics and Computer Science
Clarkson University, Potsdam, NY    13699-5815
phone:   315-268-2379       FAX:  315-268-2371
fulton@clarkson.edu
