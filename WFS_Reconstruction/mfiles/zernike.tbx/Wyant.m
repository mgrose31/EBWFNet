function [l,n]=Wyant(index)
% SYNTAX:
% [l,n]=Wyant(index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% J. C. Wyant & K. Creath,
% "Basic Wavefront Aberration Theory for Optical Metrology,
% "APPLIED OPTICS AND OPTICAL ENGINEERING, VOL. Xl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% index [ ] = 1D Zernike index (starting with tilts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% n [ ] = radial order (starting at 1 for tilts)
% l [ ] = angular order (one of -n:2:n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: Wyant.m 3036 2010-09-23 21:22:50Z amoran $

%% BEGIN_CODE 

r=index+1; % Start with tilt, not piston
t = 2*ceil(-1 + sqrt(1 + (r-1) ) ); % Total order, t = n+abs(l)
i_off = t*(t+4)/4 - index; % Number from final index with this 't'
n = t - ceil(i_off/2); % Radial order t/2 <= n <= t (in twos)
l = (t-n)*(-1)^(1-i_off); % Angular order abs(l) <= t/2
% disp(['i = ',int2str(index),' ; t = ',int2str(t),' ; ioff = ',int2str(i_off),' ; n = ',int2str(n),' ; l = ',int2str(l)])

return;