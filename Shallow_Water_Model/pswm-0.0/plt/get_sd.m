function [t,smin,smax,savg,fcor,beta,cval,lfld,labelf] = get_sd( runid, rundir )
% GET_SD  Gets scalar diagnostics from file runid.out written by pswm
%
% Usage:  [t,smin,smax,savg,fcor,beta,cval,lfld,labelf] = get_sd(runid,rundir)
%
% Input:
%     runid  run identifier [default:  see get_runid]
%     rundir directory where output files are [default:  '../run/']
%
% Output:
%     t      vector of output times
%     smin   smin(i,ifld) = minimum value of field  ifld  at time  t(i) 
%     smax   smax(i,ifld) = maximum value of field  ifld  at time  t(i) 
%     savg   savg(i,ifld) = average value of field  ifld  at time  t(i)
%     fcor   Coriolis parameter
%     beta   y-derivative of Coriolis parameter [=zero for pswm]
%     cval   phase speed
%     lfld   vector of one-character field labels
%     labelf vector of longer field labels
%
% See also:  PLOT_SD, GET_FIELD, PLOT_FIELD

% Author:  Scott R. Fulton
%
% Revision history:
% 29 Dec 2005 original version (for cswm)
% 08 Aug 2007 updated for pswm

if (~exist('rundir','var') | isempty(rundir)); rundir = '../run/'; end
if (~exist('runid' ,'var') | isempty(runid )); runid  = get_runid(rundir); end
fname = [rundir runid '.out' ];

fid = fopen(fname);

% read header line 1
model = fscanf(fid,'%s',1);
junk = fscanf(fid,'%s',1);
vernum  = fscanf(fid,'%g',1);
verdate = fscanf(fid,'%s',1);

% read header line 2
junk = fscanf(fid,'%s',1); fcor = fscanf(fid,'%g',1);
junk = fscanf(fid,'%s',1); beta = fscanf(fid,'%g',1);
junk = fscanf(fid,'%s',1); cval = fscanf(fid,'%g',1);

i = 0;
smin = [];
smax = [];
savg = [];
while 1 % loop on the time step

% end of file?
if fgetl(fid)==-1, break, end
if fgetl(fid)==-1, break, end
i = i+1;

% read the output time
junk = fscanf(fid,'%s',1); step(i) = fscanf(fid,'%d',1);
junk = fscanf(fid,'%c',8);    t(i) = fscanf(fid,'%g',1);
junk = fscanf(fid,'%s',3);

% read the scalar diagnostics for this output time
nfld  = 10;
dvals = zeros(3,nfld);
lfld  = []; labelf  = [];
for ifld=1:nfld
   junk = fscanf(fid,'%c',1);
   junk = fscanf(fid,'%1c',1);
   lfld = strvcat(lfld,junk);
   junk = fscanf(fid,'%c',3);
   junk = fscanf(fid,'%1c',26);
   labelf = strvcat(labelf,junk);
   dvals(:,ifld) = fscanf(fid,'%g',3);
end % of loop on ifld

smin = [smin; dvals(1,:)];
smax = [smax; dvals(2,:)];
savg = [savg; dvals(3,:)];

end % of loop on the time step
fclose(fid);

% turn things into column vectors
t = t';
