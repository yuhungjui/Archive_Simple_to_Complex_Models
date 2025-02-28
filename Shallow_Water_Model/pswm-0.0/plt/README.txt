File README.txt for directory $(PSWM)/plt

Periodic Shallow Water Model:  Matlab codes for plotting the model output

Files:
-----
show_run.m	shows a movie of one field from a model run
plot_run.m	makes contour plots all fields in a given model run
plot_sd.m	plots scalar diagnostics for a model run
get_sd.m	reads scalar diagnostics for a model run
plot_field.m	plots a single field
get_field.m	reads a field file written by pswm
get_cs.m	gets a cross-section of a field (at x or y value or mean)
tstring.m	constructs a time string for output

Documentation:
-------------
For each of these functions you can (from within Matlab) type "help <name>" 
to display a short description of what it does and how to use it.

Note:
----
I wrote these for testing the model and find them useful, but they may not be
particularly understandable.  And if you don't have Matlab available....
