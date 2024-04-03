function [w, s] = wpgzernw(wpg, radius, mask, nz, ns)
% SYNTAX:
% [w, s] = wpgzernw(wpg, radius, mask, nz, ns)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% wpg [ ] = 
% radius [ ] = 
% mask [ ] = 
% nz [ ] =
% ns [ ] =
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% w [ ] = 
% s [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: wpgzernw.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

s.outradius = radius;
s.spacing = min([wpg.dx, wpg.dy]);
nxy = wpg.nx * wpg.ny;
s.xsub = [];
s.ysub = [];
xmin = inf;
xmax = -inf;
ymin = inf;
ymax = -inf;
for jj=1:wpg.nx
   for ii=1:wpg.ny
      if (mask(jj,ii) > 0)
         xmin = min([xmin,wpg.u(jj)]);
         xmax = max([xmax,wpg.d(jj)]);         
         s.xsub(end+1) = wpg.u(jj);
         ymin = min([ymin,wpg.l(ii)]);
         ymax = max([ymax,wpg.r(ii)]);
         s.ysub(end+1) = wpg.l(ii);
      end;
   end;
end;
s.xcenter = (xmax + xmin)/2;
s.ycenter = (ymax + ymin)/2;
w = mkwv(s,nz,ns);
