function [f,x,y,t,fcor,beta,cval,lab,runid] = get_field( fname )
% GETFLD  Gets a field from a file written by pswm
%
% Usage:  [f,x,y,t,fcor,beta,cval,lab,runid] = get_field( fname )
%
% Input:
%     fname  must be a string constant, e.g., get_field( 'test_z.t01' )
%
% Output:
%     f     field read [the remaining variables are optional]
%     x     corresponding x values
%     y     corresponding y values
%     t     corresponding time (scalar)
%     fcor  Coriolis parameter
%     beta  y-derivative of Coriolis parameter
%     cval  phase speed
%     lab   one-character field label (one of:  u, v, p, d, z) 
%     runid string identifying the run
%
% See also:  PLOT_RUN, PLOT_FIELD, GET_SD, PLOT_SD

% Author:  Scott R. Fulton
%
% Revision history:
% 09 Dec 2005 original version (for cswm)
% 29 Dec 2005 minor updates to match changes in cswm

fid = fopen(fname);

% read header line 1
junk = fscanf(fid,'%s',2);
vernum  = fscanf(fid,'%g',1);
verdate = fscanf(fid,'%s',1);

% read header line 2
junk = fscanf(fid,'%s',1);
lab  = fscanf(fid,'%s',1);
junk = fscanf(fid,'%s',2);
runid = fscanf(fid,'%s',1);

% read header line 3
nxny = fscanf(fid,'%g',2);
nx = nxny(1);
ny = nxny(2);

% read header line 4
lims = fscanf(fid,'%g',4);
xa = lims(1);
xb = lims(2);
ya = lims(3);
yb = lims(4);

% read header line 5
vals = fscanf(fid,'%g',4);
fcor = vals(1);
beta = vals(2);
cval = vals(3);
t    = vals(4);

% create the coordinate arrays
x = linspace(xa,xb,nx+1);
y = linspace(ya,yb,ny+1);

% read the field values
f = zeros(nx+1,ny+1);
f = fscanf(fid,'%g',[nx+1,ny+1]);
fclose(fid);
