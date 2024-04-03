function wg = wfspixgeom2d(nx, x0, dx, pxb, pxf, xmax, ny, y0, dy, pyb, ...
              pyf, ymax, intlevel);
% SYNTAX: 
% wg = wfspixgeom2d(nx, x0, dx, pxb, pxf, xmax, ny, y0, dy, pyb, pyf, ymax, 
%      intlevel);
% wg = wfspixgeom2d(nx, xy0, dxy, pxyb, pxyf, xymax, intlevel);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: wfspixgeom2d.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

if ((nargin == 6) | (nargin == 12))
   intlevel = [];
elseif (nargin == 7)
   intlevel = ny;
elseif (nargin ~= 13)
   return;
end;
wgx = wfspixgeom1d(nx, x0, dx, pxb, pxf, xmax, intlevel);
if (nargin >= 12)
   wgy = wfspixgeom1d(ny, y0, dy, pyb, pyf, ymax, intlevel);
else
   wgy = wgx;
end;
wg = wgx;
wg.y = wgy.x;
wg.dy = wgy.dx;
wg.l = wgy.u;
wg.pyb = wgy.pxb;
wg.r = wgy.d;
wg.pyf = wgy.pxf;
wg.npy = wgy.npx;
wg.ny = wgy.nx;
wg.ymax = wgy.xmax;
wg.yintlevel = wgy.xintlevel;

