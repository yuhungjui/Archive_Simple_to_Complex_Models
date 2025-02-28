function plot_run( fields, runid, rundir, times, nrows )
% PLOT_RUN  Plots the solution from a pswm run
%
% Usage:  plot_run( fields, runid, rundir, times, nrows )
%
% Input:
%     fields list of fields to plot (see possible values in plot_field)
%            example:  fields = ['W'; 'z']
%            default:  plot every variable available
%     runid  run identifier [default:  see get_runid]
%     rundir directory where output files are [default:  ../run]
%     times  list of times to plot (indices of files) [default:  all]
%     nrows  number of rows to plot per page [default:  matches columns]
%
% Output:
%     Plots specified fields:  one time per row, one field per column
%
% See also:  PLOT_FIELD, GET_FIELD

% Author:  Scott R. Fulton
%
% Revision history:
% 10 Dec 2005 original version (for cswm)
% 06 Jan 2006 changed order of arguments
% 08 Aug 2007 updated for pswm

if (~exist('rundir','var') | isempty(rundir) ), rundir = '../run/'; end
if (~exist('runid' ,'var') | isempty(runid ) ), runid  = get_runid(rundir); end

% determine which variables are available
if (~exist('fields','var') | isempty(fields) )
fname = [rundir, runid '_u.t00' ]; plotu = (exist(fname,'file')==2);
fname = [rundir, runid '_v.t00' ]; plotv = (exist(fname,'file')==2);
fname = [rundir, runid '_p.t00' ]; plotp = (exist(fname,'file')==2);
fname = [rundir, runid '_d.t00' ]; plotd = (exist(fname,'file')==2);
fname = [rundir, runid '_z.t00' ]; plotz = (exist(fname,'file')==2);
plotw =  (plotu & plotv); if plotw, plotu=0; plotv=0; end
plotW =  (plotw & plotp); if plotW, plotw=0; plotp=0; end
fields = [];
if plotW; fields = [fields; 'W']; end
if plotw; fields = [fields; 'w']; end
if plotu; fields = [fields; 'u']; end
if plotv; fields = [fields; 'v']; end
if plotp; fields = [fields; 'p']; end
if plotd; fields = [fields; 'd']; end
if plotz; fields = [fields; 'z']; end
end

if (~exist('times','var')), times = []; end

% determine the number of rows and columns to plot
ncols = length(fields);
if ncols==0; error(['no variables to plot for runid=' runid, ...
                    ' in rundir=' rundir]); end

if (~exist('nrows','var') | isempty(nrows) ), nrows = ncols; end
if nrows<1, nrows = 1; end

% loop over output time
clf
itime = -1;
irow  = 0; iplot = 0;
while true
   itime = itime+1;
   % quit if no more files are present 
   if (length(dir([rundir runid '_*.t' num2str(itime,'%02d')]))==0), break; end
   % skip this time if it wasn't specified
   if (length(times)>0 & ~length(find(times==itime))), continue; end
   irow = irow+1;
   if irow>nrows % time to start a new page
      pause; clf; irow = 1; iplot = 0;
   end
   labels = [0 2 3]; if irow==nrows, labels(1) = 2; end
   for ifld=1:length(fields)
      iplot = iplot+1; subplot(nrows,ncols,iplot);
      plot_field( fields(ifld), runid, rundir, itime, labels);
      labels(2) = 0;
   end
end % of loop over output time
