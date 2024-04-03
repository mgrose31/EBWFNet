function [s, p] = reconcloop(gam, gaminv, os, n, nf);
% SYNTAX:
% [s, p] = reconcloop(gam, gaminv, os, n, nf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% gam [ ] = 
% gaminv [ ] = 
% os [ ] = 
% n [ ] = 
% nf [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% s [ ] = 
% p [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: reconcloop.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

ns = prod(size(os));
s = zeros(size(gam,1),n);
p = zeros(size(gaminv,1),n);
snow = os + nf*mean(abs(os))*randn(size(os));
for i=1:n
   p(:,i) = gaminv*snow;
   s(:,i) = gam*p(:,i);
   snow = snow - s(1:ns,i);
   if (nf ~= 0); snow = snow + nf*mean(abs(snow))*randn(size(snow)); end;
end;