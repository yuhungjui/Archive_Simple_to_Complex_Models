function [fcs,xy,t,fcor,beta,cval,lab,runid] = get_cs( fname, xy, valxy )
% GET_CS  Gets a cross-section of a field from a file written by pswm
%
% Usage:  [fcs,xy,t,fcor,beta,cval,lab,runid] = get_cs( fname, xy, valxy )
%
% Input:
%     fname  must be a string constant, e.g., get_cs( 'test_z.t01' )
%     xy     specifies the cross-section:
%            xy='x':  cross-section in x at y=valxy (or mean in y)
%            xy='y':  cross-section in y at x=valxy (or mean in x) [default]
%     valxy  x or y value for cross-section (or 'mean' for mean [default])
%
% Output:
%     fcs   cross-section of field [the remaining variables are optional]
%     xy    corresponding x or y values
%     t     corresponding time (scalar)
%     fcor  Coriolis parameter
%     beta  y-derivative of Coriolis parameter
%     cval  phase speed
%     lab   one-character field label (one of:  u, v, p, d, z) 
%     runid string identifying the run
%
% See also:  PLOT_CS, GET_FIELD

% Author:  Scott R. Fulton
%
% Revision history:
% 31 Dec 2005 original version (for cswm)
% 08 Aug 2007 updated for pswm

% read the field
[f,x,y,t,fcor,beta,cval,lab,runid] = get_field( fname );

% construct the cross-section
if (~exist('xy' ,'var') | isempty(xy)); xy = 'x';  end
if (~exist('valxy' ,'var') | isempty(valxy)); valxy = 'mean';  end
switch lower(xy)
   case 'x' % cross-section in x at y=valxy (or mean in y)
      xy = x';
      if valxy=='mean'
         fcs = mean(f,2);
      else 
         % find the interval containing the specified value
         dy = y(2)-y(1);
         j  = find(abs(y-valxy)<abs(dy));
         if length(j)==1 % copy values at this point
            fcs = f(:,j);
         elseif length(j)==2 % linear interpolation
            fcs = ((y(j(2))-valxy)*f(:,j(1)) + (valxy-y(j(1)))*f(:,j(2)))/dy;
         else
            y
            valxy
            error('cannot find y=valxy');
         end
      end
   otherwise % cross-section in y at x=valxy (or mean in x)
      xy = y';
      if valxy=='mean'
         fcs = mean(f,1);
      else 
         % find the interval containing the specified value
         dx = (2)-(1);
         i  = find(abs(x-valxy)<abs(dx));
         if length(i)==1 % copy values at this point
            fcs = f(i,:);
         elseif length(i)==2 % linear interpolation
            fcs = ((x(i(2))-valxy)*f(i(1),:) + (valxy-x(i(1)))*f(i(2),:))/dx;
         else
            x
            valxy
            error('cannot find x=valxy');
         end
         fcs = fcs'; % to get a column vector
      end
end % of switch statement for x or y cross-section
