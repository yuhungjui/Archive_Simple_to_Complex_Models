function time = tstring( t, first, last )
% TSTRING  Returns string expressing time t in the form dd:hh:mm:ss
%
% Usage:  time = tstring( t, first, last )
%
% Input:
%     t      time in seconds
%     first  [optional] first unit:  one of 'd', 'h', 'm', 's'		
%     last   [optional] last  unit:  one of 'd', 'h', 'm', 's'		
%
% Output:
%     time   character string of the form  dd:hh:mm:ss
%
% Examples:  
%     tstring( 117 )           returns '00:00:01:57'
%     tstring( 117, 'm' )      returns '01:57'
%     tstring( 117, 'h', 'm' ) returns '00:02'

% Author:  Scott R. Fulton
%
% Revision history:
% 12/10/05 original version

if (~exist('first','var') | isempty(first) ), first = 'd'; end
if (~exist('last' ,'var') | isempty(last ) ), last  = 's'; end

% determine the number of seconds, minutes, hours, and days
%s = round(t)
%m = round((t-s)/60)
%h = round((t-s-m*60)/(60*60))
%d = round((t-s-m*60-h*60*60)/(24*60*60))
minute = 60; hour = 60*minute; day = 24*hour;
d = floor(t/day);
h = floor((t-d*day)/hour);
m = floor((t-d*day-h*hour)/minute);
s = t-d*day-h*hour-m*minute;

% by default, diplay all units
sd = 1; sh = 1; sm = 1; ss = 1;

% correct for the last value to be displayed
switch lower(last)
   case {'m'}
      m = round(m+s/minute);
      if m==60, m=0; h=h+1; end;
      if h==24, h=0; d=d+1; end;
      ss = 0;
   case {'h'}
      h = round(h+(m*minute+s)/hour); 
      if h==24, h=0; d=d+1; end;
      sm = 0; ss = 0;
   case {'d'}
      d = round(d+(h*hour+m+minute+s)/day); 
      sh = 0; sm = 0; ss = 0;
   otherwise
      s = round(s);
      if s==60, s=0; m=m+1; end;
      if m==60, m=0; h=h+1; end;
      if h==24, h=0; d=d+1; end;
end

% correct for the first value to be displayed
switch lower(first)
   case {'h'}
      h = h + d*24;
      sd = 0;
   case {'m'}
      m = m + (h + d*24)*60;
      sd = 0; sh = 0;
   case {'s'}
      s = s + (m + (h + d*24)*60)*60;
      sd = 0; sh = 0; sm = 0;
end

% construct the string
time = [];
if sd
   time = [time,sprintf('%02d',d)];
end
if sh
   if (~isempty(time)), time = [time,':']; end
   time = [time, sprintf('%02d',h)]; 
end
if sm
   if (~isempty(time)), time = [time,':']; end
   time = [time, sprintf('%02d',m)]; 
end
if ss
   if (~isempty(time)), time = [time,':']; end
   time = [time, sprintf('%02d',s)]; 
end
