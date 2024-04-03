function [dm, x, y] = actmap(d, xact, yact, bval, reltol, padact)
% SYNTAX:
% [dm, x, y] = actmap(d, xact, yact, bval, reltol, padact)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% d [ ] = 
% xact [ ] = 
% yact [ ] = 
% bval [ ] = 
% reltol [ ] = 
% padact [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% dm [ ] = 
% x [ ] = 
% y [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: actmap.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

if (~exist('bval','var'))
   bval = nan;
end;
if (isempty(bval))
   bval = nan;
end;
if (~exist('reltol','var'))
   reltol = 1.0e-06;
end;
if (isempty(reltol))
   reltol = 1.0e-06;
end;
nc = prod(size(xact));
dxs = abs(xact(2:nc)-xact(1:(nc-1)));
rdxs = reltol * max(dxs);
idxs = find(dxs > rdxs);
dxs = dxs(idxs);
dx = min(dxs);
dys = abs(yact(2:nc)-yact(1:(nc-1)));
rdys = reltol * max(dys);
idys = find(dys > rdys);
dys = dys(idys);
dy = min(dys);
eps = reltol * min([dx,dy]);
x = xact(1);
y = yact(1);
for i=2:nc
   if (length(find(abs(x-xact(i)) < eps)) <= 0)
      x = [x, xact(i)];
   end;
   if (length(find(abs(y-yact(i)) < eps)) <= 0)
      y = [y, yact(i)];
   end;
end;
x = sort(x);
y = sort(y);
nx = prod(size(x));
ny = prod(size(y));
dm = repmat(bval,nx,ny);
for i=1:nc
   ix = find(abs(x-xact(i)) < eps);
   iy = find(abs(y-yact(i)) < eps);
   dm(ix,iy) = d(i);
end;
if (~exist('padact','var'))
   padact = 0;
end;
if (isempty(padact))
   padact = 0;
end;
for ii=1:padact
   x = [x(1)-(x(2)-x(1)), x, x(end)+(x(end)-x(end-1))];
   y = [y(1)-(y(2)-y(1)), y, y(end)+(y(end)-y(end-1))];
end;
if (padact > 0)
   dmhold = dm;
   dm = repmat(nan, nx+2*padact, ny+2*padact);
   dm(1+padact:end-padact,1+padact:end-padact) = dmhold;
end;
return
