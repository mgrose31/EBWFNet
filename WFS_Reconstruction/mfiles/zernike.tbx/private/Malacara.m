function [l,n]=Malacara(index)
% SYNTAX:
% [l,n]=Malacara(index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: Malacara.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE 

%function [l,n]=Malacara(index)
% Malacara 1D ordering
% Inputs:
% index - 1D Zernike index (starting with tilts)
%
% Outputs:
% n - radial order (starting at 1 for tilts)
% l - angular order (one of -n:2:n)

r=index+1;
n=ceil((-3 + sqrt(9 + 8*(r-1)))/2);
l = n^2 + 2*(n - r + 1);

return;