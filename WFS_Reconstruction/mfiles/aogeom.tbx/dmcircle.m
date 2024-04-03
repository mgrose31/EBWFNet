function handle = dmcircle(coords, x, y, color)
% SYNTAX:
% handle = dmcircle(coords, x, y, color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% coords [ ] = 
% x [ ] = 
% y [ ] = 
% color [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% handle [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: dmcircle.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE

handle = zeros(2);
outerH = d_circle(coords(2), x, y,color);
if(coords(1) > 0)
   innerH = d_circle(coords(1), x, y, color);
else
   innerH = -1;
end	
handle(1) = innerH;
handle(2) = outerH;
