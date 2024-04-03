function [pm, p, x, y] = wafflep(c, tol)
% SYNTAX:
% [pm, p, x, y] = wafflep(c, tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%y
% INPUTS:
% c [ ] =
% tol [ ] =
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% pm [ ] = 
% p [ ] = 
% x [ ] = 
% y [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: wafflep.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

if (nargin < 2)
   tol = 0.001;
end;
if (~exist('tol','var'))
   tol = 0.001;
end;
if (isempty(tol))
   tol = 0.001;
end;   
sc = size(c.x);
nc = prod(sc);
x = c.x(1);
y = c.y(1);
for i=2:nc
   if (length(find(abs(x-c.x(i)) < tol)) <= 0)
      x = [x, c.x(i)];
   end;
   if (length(find(abs(y-c.y(i)) < tol)) <= 0)
      y = [y, c.y(i)];
   end;
end;
x = sort(x);
y = sort(y);
nx = length(x);
ny = length(y);
p = zeros(nc,1);
pm = zeros(nx,ny);
for i=1:nc
   ix = find(abs(x-c.x(i)) < tol);
   iy = find(abs(y-c.y(i)) < tol);
   jx = -1;
   if (mod(ix,2) == 0); jx = 1; end;
   jy = -1;
   if (mod(iy,2) == 0); jy = 1; end;
   p(i) = jx*jy;
   pm(ix,iy) = p(i);
end;
return
