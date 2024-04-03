function [ic] = findcorner(c, x, y, tolin);
% SYNTAX:
% [ic] = findcorner(c, x, y, {tolin})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: findcorner.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

if (nargin > 3)
   tol = tolin;
else
   relsize = min([x,y]);
   if (x == 0) relsize = y; end;
   if (y == 0) relsize = x; end;
   if (relsize == 0) relsize = 0.1; end;
   tol = relsize * 0.01;
end;
sc = size(c.ns);
nc = prod(sc);
ic = 1;
while ((ic <= nc) & (c.ns(ic) > 0))
   if ((abs(c.x(ic)-x) < tol) & (abs(c.y(ic)-y) < tol))
      return;
   end;
   ic = ic + 1;
end;
if (ic > nc) ic = -1;
return;
end
