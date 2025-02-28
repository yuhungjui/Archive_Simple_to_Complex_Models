function plot_field( field, runid, rundir, nfile, labels, ptype )
% PLOT_FIELD  Plots a single field from a pswm output file
%
% Usage:  plot_field( field, runid, rundir, nfile, labels, ptype )
%
% Input:
%     field  one-character label for field(s) to plot:
%            u  x-component of wind
%            v  y-component of wind
%            p  scaled deviation geopotential
%            d  divergence
%            z  relative vorticity
%            e  absolute  vorticity
%            q  potential vorticity
%            w  wind (vectors)
%            W  wind (vectors) and geopotential [default]
%     runid  run identifier [default:  see get_runid]
%     rundir directory where output files are [default:  '../run/']
%     nfile  index of output file [default:  0]
%     labels vector to control drawing of labels:
%            labels(1) is for x axis: 0=none, 1=scale, 2=label
%            labels(2) is for y axis: 0=none, 1=scale, 2=label
%            labels(3) is for  title: 0=none, 1=field, 2=time, 3=both
%            Default:  xylab = [2 2 3]
%     ptype  if present, specifies the plot type:
%            'fc'   filled contours [default]
%            'c'    unfilled contours
%            'mesh' 3d mesh plot
%            'surf' 3d surface plot
%            Note:  ptype is ignored for field='w' and field='W'
%
% Output:
%     Draws a contour plot of a single field (or wind vectors)
%
% Example:  
%     plot_field('z',[],'test',2) reads file  ../run/test_z.t02 and plots it
%
% See also:  PLOT_RUN, GET_FIELD

% Author:  Scott R. Fulton
%
% Revision history:
% 10 Dec 2005 original version (for cswm)
% 06 Jan 2006 changed order of arguments
% 08 Jan 2006 added 3d plots to plot type
% 08 Aug 2007 updated for pswm

if (~exist('field' ,'var') | isempty(field )); field  = 'W'; end
if (~exist('rundir','var') | isempty(rundir)); rundir = '../run/'; end
if (~exist('runid' ,'var') | isempty(runid )); runid  = get_runid(rundir); end
if (~exist('nfile' ,'var') | isempty(nfile )); nfile  = 0; end
if (~exist('labels','var') | isempty(labels)); labels = [2 2 3]; end
if (~exist('ptype' ,'var') | isempty(ptype )); ptype = 'fc'; end

% read the input file(s) and construct the field(s)
switch field
   case {'u','v','p','d','z'} % primary variables output by pswm
       fname = [rundir runid '_' field '.t' sprintf('%02d',nfile) ];
       [f,x,y,t,fcor,beta,cval] = get_field( fname );
   case {'e','q'} % absolute and potential vorticity
       fname = [rundir runid '_z.t' sprintf('%02d',nfile) ];
       [f,x,y,t,fcor,beta,cval] = get_field( fname );
       for j=1:length(y)
          f(:,j) = (fcor + beta*y(j)) + f(:,j);
       end
       if field=='q'
          fname = [rundir runid '_p.t' sprintf('%02d',nfile) ];
          [p,x,y,t,fcor,beta,cval] = get_field( fname );
          f = f./(1 + p/cval);
       end
   case {'w','W'} % wind vectors
       fname = [rundir runid '_u.t' sprintf('%02d',nfile) ];
       [u,x,y,t,fcor,beta,cval] = get_field( fname );
       fname = [rundir runid '_v.t' sprintf('%02d',nfile) ];
       [v,x,y,t,fcor,beta,cval] = get_field( fname );
       if field=='W'
          fname = [rundir runid '_p.t' sprintf('%02d',nfile) ];
          [p,x,y,t,fcor,beta,cval] = get_field( fname );
       end
   otherwise
      error(['Unknown field:  ',field]);
end

% scale the coordinates
x = x/1000; y = y/1000;

% draw the plot
switch field
   case {'w'} % wind vectors
      quiver(x,y,u',v'); axis manual;
      axis image
      axis([min(x) max(x) min(y) max(y)]);
      lfield = 'wind';
   case {'W'} % wind vectors and geopotential
      %contourf(x,y,p',cval);
      %cval = [-20:2:20]; contourf(x,y,p',cval); 
      contourf(x,y,p'); 
      hold on; axis image
      quiver(x,y,u',v'); hold off;
      lfield = 'wind and p';
   otherwise
      switch ptype
         case {'c'}
            contour(x,y,f');
            colorbar
            axis image
         case {'mesh'}
            mesh(x,y,f');
         case {'surf'}
            surf(x,y,f');
         otherwise
            contourf(x,y,f');
            colorbar
            axis image
      end
      lfield = field;
end

% draw the labels
if (labels(1)<1), set(gca,'XTick',[]); end
if (labels(2)<1), set(gca,'YTick',[]); end
if (labels(1)>1), xlabel('x (km)'); end
if (labels(2)>1), ylabel('y (km)'); end
switch labels(3)
   case 1
      title([lfield ' field']); 
   case 2
      title(['t=' tstring(t,'h')]); 
   case 3
      title([lfield ' at t=' tstring(t,'h')]); 
end
