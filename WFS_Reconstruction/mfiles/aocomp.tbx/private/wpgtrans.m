function wg = wpgtrans(wgin, movedown, moveright)
% SYNTAX: 
% wg = wpgtrans(wgin, movedown, moveright)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: wpgtrans.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

wg = wgin;
wg.x = wg.x + movedown;
if (wg.xintlevel > 1)
   wg.x = round(wg.x);
end;
ind = find(wg.x > wg.xmax);
if (~isempty(ind))
   wg.x(ind) = wg.xmax;
end;
%
wg.dx = wg.x(2:wg.nx) - wg.x(1:(wg.nx-1));
%
wg.u = wg.x - wg.pxb;
if (wg.xintlevel > 0)
   wg.u = round(wg.u);
end;
ind = find(wg.u < 1);
if (~isempty(ind))
   wg.u(ind) = 1;
end;
wg.pxb = wg.x - wg.u;
%
wg.d = wg.x + wg.pxf;
if (wg.xintlevel > 0)
   wg.d = round(wg.d);
end;
ind = find(wg.d > wg.xmax);
if (~isempty(ind))
   wg.d(ind) = wg.xmax;
end;
wg.pxf = wg.d - wg.x;
%
wg.npx = wg.d - wg.u + 1;
%
if (~exist('moveright','var') | ~isfield(wg, 'y'))
   return;
end;
%
wg.y = wg.y + moveright;
if (wg.yintlevel > 1)
   wg.y = round(wg.y);
end;
ind = find(wg.y > wg.ymax);
if (~isempty(ind))
   wg.y(ind) = wg.ymax;
end;
%
wg.dy = wg.y(2:wg.ny) - wg.y(1:(wg.ny-1));
%
wg.l = wg.y - wg.pyb;
if (wg.yintlevel > 0)
   wg.l = round(wg.l);
end;
ind = find(wg.l < 1);
if (~isempty(ind))
   wg.l(ind) = 1;
end;
wg.pyb = wg.y - wg.l;
%
wg.r = wg.y + wg.pyf;
if (wg.yintlevel > 0)
   wg.r = round(wg.r);
end;
ind = find(wg.r > wg.ymax);
if (~isempty(ind))
   wg.r(ind) = wg.ymax;
end;
wg.pyf = wg.r - wg.y;
%
wg.npy = wg.r - wg.l + 1;
