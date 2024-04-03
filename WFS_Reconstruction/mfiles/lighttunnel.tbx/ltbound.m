function [lx, rx, minx, maxx, cxx] = ltbound(cx, tx, bx);
% SYNTAX:
% [lx, rx, minx, maxx, cxx] = ltbound(cx, tx, bx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Determines width of add-on 'buffer' to image for lighttunnel.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% cx [ ] = spatial coordinates of image
% tx [ ] = tilt function per pixel (same size as image)
% bx [ ] = blur function per pixel (same size as image)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% lx [ ] = 'left' border function
% rx [ ] = 'right' border function
% minx [ ] = lowest value of coordinates for image
% maxx [ ] = highest value of coordinates for image
% cxx [ ] = coordinate array including 3rd dimension of input image array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEE ALSO: lighttunnel.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Robert W. Praus, II
% (c) 2005 MZA Associates Corporation, Albuquerque, NM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: ltbound.m 3036 2010-09-23 21:22:50Z amoran $

%% Tilt and Weighted Blurring Boundaries
lx = tx - bx *0.5;                   %   lx (left) =tilt - weighted blurring per pixel
rx = tx + bx *0.5;                  %   rx (right) =tilt + weighted blurring per pixel
minx = floor(min(lx(:)));%  minimum tilt - weighted blurring factor per pixel
maxx = ceil(max(rx(:)));%   maximum tilt + weighted blurring factor per pixel

%% Make 3-d set of spatial coordinates per pixel
cxx = repmat(cx,[1,1,size(lx,3)]);
lx = cxx + lx;
rx = cxx + rx;
