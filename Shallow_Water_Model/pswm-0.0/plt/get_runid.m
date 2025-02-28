function runid = get_runid( rundir )
% GET_RUNID  Gets the runid from the last run in rundir
%
% Usage:  runid = get_runid( rundir )
%
% Input:
%     rundir run directory [default:  '../run/']
%
% Output:
%     runid  run identifier of the last run (from the file runid.txt)
%
% See also:  PLOT_RUN, PLOT_FIELD, PLOT_SD, ...

% Author:  Scott R. Fulton
%
% Revision history:
% 02 Jan 2006 original version (for cswm)
% 15 Aug 2007 updated for pswm (fixed the comments)

if (~exist('rundir','var') | isempty(rundir)); rundir = '../run/'; end

% read the run identifier from the file runid.txt
fname = [rundir 'runid.txt' ];
fid = fopen(fname);
runid = fscanf(fid,'%s',1);
fclose(fid);
