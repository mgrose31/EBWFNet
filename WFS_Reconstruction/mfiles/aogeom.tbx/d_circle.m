function handle = d_circle(radius, offset_x, offset_y, color)
% SYNTAX:
% handle = d_circle(radius, offset_x, offset_y, color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% radius [ ] = 
% offset_x [ ] =
% offset_y [ ] = 
% color [ ] =
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% handle [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: d_circle.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE

numpts = 60;
t = 0:(2*pi/numpts):(2*pi);
x = cos(t)*radius + offset_x;
y = sin(t)*radius + offset_y;
handle = line(x,y,'color',color);

