function [pm, x, y] = lsfmap(p, c, bval, reltol)
% SYNTAX:
% [pm, x, y] = lsfmap(p, c, bval, reltol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% p [ ] = 
% c [ ] = 
% bval [ ] = 
% reltol [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% pm [ ] = 
% x [ ] = 
% y [ ] =
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: lsfmap.m 3061 2010-10-07 21:13:39Z amoran $

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
bval = min(min(bval));
sc = size(c.ns);
nc = prod(sc);
x = sort(unique(c.x));
y = sort(unique(c.y));
nx = prod(size(x));
ny = prod(size(y));
rdx = reltol * max(abs(x(2:nx)-x(1:nx-1)));
rdy = reltol * max(abs(y(2:ny)-y(1:ny-1)));
[ws,wf]=warning; warning off;
dx = abs((x(2:nx) - x(1:nx-1))./x(1:nx-1));
warning(ws); warning(wf);
ind = find(dx > rdx);
ni = prod(size(ind));
if (ni ~= nx)
   x = [x(1),x(ind+1)];
   nx = prod(size(x));
end;
[ws,wf]=warning; warning off;
dy = abs((y(2:ny) - y(1:ny-1))./y(1:ny-1));
warning(ws); warning(wf);
ind = find(dy > rdy);
ni = prod(size(ind));
if (ni ~= ny)
   y = [y(1),y(ind+1)];
   ny = prod(size(y));
end;
pm = repmat(bval,nx,ny);
for i=1:nc
   if (abs(c.x(i)) < rdx)
      ix = find(abs(x) < rdx);
   else
      ix = find(abs((x-c.x(i))./c.x(i)) < rdx);
   end;
   if (abs(c.y(i)) < rdy)
      iy = find(abs(y) < rdy);
   else
      iy = find(abs((y-c.y(i))./c.y(i)) < rdy);
   end;
   pm(ix,iy) = p(i);
end;
return
