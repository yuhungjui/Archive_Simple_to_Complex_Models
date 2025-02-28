function M = show_run( field, runid, rundir, flim, times )
% SHOW_RUN  Shows a movie of the solution from a pswm run
%
% Usage:  M = show_run( field, runid, rundir, flim, times )
%
% Input:
%     field  field to plot (see possible values in pltfld) [default:  'z']
%     runid  run identifier [default:  see get_runid]
%     rundir directory where output files are [default:  ../run]
%     flim   scale limits for field [default:  [-Inf Inf] for autoscale]
%     times  list of times to plot (indices of files) [default:  all]
%
% Output:
%     M     Movie file to be played with movie(M)
%           (show_run plays this once as it constructs the movie)
%
% See also:  PLOT_RUN

% Author:  Scott R. Fulton
%
% Revision history:
% 31 Dec 2005 original version (for cswm)
% 02 Jan 2006 fixed axes so x and y are properly represented
% 06 Jan 2006 changed order of arguments
% 08 Aug 2007 updated for pswm

if (~exist('rundir','var') | isempty(rundir) ), rundir = '../run/'; end
if (~exist('runid' ,'var') | isempty(runid ) ), runid  = get_runid(rundir); end

% determine which field and times to plot
if (~exist('field','var') | isempty(field) ), field = 'z'; end
if (~exist('times','var')), times = []; end

% get the field data
itime = -1; j = 0;
while true % loop over output time
   itime = itime+1;
   % quit if no more files are present 
   if (length(dir([rundir runid '_*.t' num2str(itime,'%02d')]))==0), break; end
   % skip this time if it wasn't specified
   if (length(times)>0 & ~length(find(times==itime))), continue; end
   fname = [rundir runid '_' field '.t' sprintf('%02d',itime) ];
   [f,x,y,t,fcor,beta,cval] = get_field( fname );
   switch field
      case {'u','v','p'}
         f = f/cval;
         labelf = [field '/c'];
      case {'d','z'}
         if fcor>0, f = f/fcor; end;
         labelf = [field '/f'];
   end
   j = j+1;
   F(:,:,j) = f;
end % of loop over output time

% determine the scale limits
x = x/1000; y = y/1000;
if (~exist('flim','var') | isempty(flim) )
   Fmin = min(min(min(F)));
   Fmax = max(max(max(F)));
   flim = [Fmin Fmax]; 
end

% make the movie
clf
nframes = j;
for j=1:nframes
   mesh(x,y,F(:,:,j)');
   xlabel('x  (km)');
   ylabel('y  (km)');
   zlabel(labelf);
   axis([-Inf Inf -Inf Inf flim]);
   M(j) = getframe;
end % of loop over frame
