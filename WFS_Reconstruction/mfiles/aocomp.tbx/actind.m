function [ix] = actind(xact, xopd, releps)
% SYNTAX:
% [ix] = actind(xact, xopd, releps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% xact [ ] = 
% xopd [ ] = 
% releps [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% ix [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: actind.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

sx = size(xact);
n = prod(sx);
dxs = abs(xact(2:n)-xact(1:(n-1)));
idxs = find(dxs > 0);
dxs = dxs(idxs);
dx = min(dxs);
if (~exist('releps','var'))
   releps = 1.0e-06;
end;
if (isempty(releps))
   releps = 1.0e-06;
end;
eps = releps * dx;
ix = zeros(sx);
for ii=1:n
   k = find(abs(xact(ii)-xopd) < eps);
   if (length(k) == 0)
      ix = ii;
      return;
   end;
   ix(ii) = k;
end;
