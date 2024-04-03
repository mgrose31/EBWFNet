function wg = wfspixgeom1d(nx, x0, dx, pxb, pxf, xmax, intlevel);
% SYNTAX: 
% wg = wfspixgeom1d(nx, x0, dx, pxb, pxf, xmax, intlevel);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: wfspixgeom1d.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

if (~exist('intlevel','var'))
   intlevel = 1; % Boundaries must be integers.
end;
if (isempty(intlevel))
   intlevel = 1; % Boundaries must be integers.
end;
%
if (prod(size(dx)) == 1)
   wg.x = [x0:dx:x0+(nx-1)*dx];
elseif (prod(size(dx)) == (nx-1))
   wg.x(1) = x0;
   for ii=2:nx
      wg.x(ii) = wg.x(ii-1) + dx(ii-1);
   end;
else
   return;
end;
if (intlevel > 1)
   wg.x = round(wg.x);
end;
ind = find(wg.x > xmax);
if (~isempty(ind))
   wg.x(ind) = xmax;
end;
wg.dx = wg.x(2:nx) - wg.x(1:(nx-1));
%
if (prod(size(pxb)) == 1)
   wg.u = wg.x - pxb;
   wg.pxb = repmat(pxb, size(wg.x));
elseif (prod(size(pxb)) == nx)
   wg.u = wg.x - pxb;
   wg.pxb = pxb;
else
   return;
end;
if (intlevel > 0)
   wg.u = round(wg.u);
end;
ind = find(wg.u < 1);
if (~isempty(ind))
   wg.u(ind) = 1;
end;
%
if (prod(size(pxf)) == 1)
   wg.d = wg.x + pxf;
   wg.pxf = repmat(pxf, size(wg.x));
elseif (prod(size(pxf)) == nx)
   wg.d = wg.x + pxf;
   wg.pxf = pxf;
else
   return;
end;
if (intlevel > 0)
   wg.d = round(wg.d);
end;
ind = find(wg.d > xmax);
if (~isempty(ind))
   wg.d(ind) = xmax;
end;
%
wg.npx = wg.d - wg.u + 1;
wg.nx = nx;
wg.xmax = xmax;
wg.xintlevel = intlevel;