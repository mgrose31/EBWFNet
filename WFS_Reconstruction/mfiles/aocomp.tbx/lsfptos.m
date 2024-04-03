function [ptos, c] = lsfptos(xsub, ysub, xact, yact)
% SYNTAX: 
% [ptos, c] = lsfptos(xsub, ysub, {xact, yact})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS: 
% xsub [ ] = 
% ysub [ ] = 
% xact [ ] = 
% yact [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% ptos [ ] = 
% c [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: lsfptos.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

ns = length(xsub);
dsubs = abs(xsub(2:ns) - xsub(1:(ns-1)));
indnz = find(dsubs);
dsub = min(dsubs(indnz));
dsubhx = dsub / 2.0;
dsubs = abs(ysub(2:ns) - ysub(1:(ns-1)));
indnz = find(dsubs);
dsub = min(dsubs(indnz));
dsubhy = dsub / 2.0;
tol = 2 * 0.1 * min([dsubhx,dsubhy]);
c.x = zeros(1,5*2*ns);
c.y = c.x;
c.ns = c.x;
c.ll = c.x;
c.lr = c.x;
c.ul = c.x;
c.ur = c.x;
maxic = 0;
for i=1:ns
   ic = findcorner(c, xsub(i)-dsubhx, ysub(i)-dsubhy, tol);
   if (c.ns(ic) == 0)
      c.x(ic) = xsub(i) - dsubhx;
      c.y(ic) = ysub(i) - dsubhy;
   end;
   c.ns(ic) = c.ns(ic) + 1;
   c.ll(ic) = i;
   maxic = max([maxic, ic]);
   ic = findcorner(c, xsub(i)+dsubhx, ysub(i)-dsubhy, tol);
   if (c.ns(ic) == 0)
      c.x(ic) = xsub(i) + dsubhx;
      c.y(ic) = ysub(i) - dsubhy;
   end;
   c.ns(ic) = c.ns(ic) + 1;
   c.lr(ic) = i;
   maxic = max([maxic, ic]);
   ic = findcorner(c, xsub(i)-dsubhx, ysub(i)+dsubhy, tol);
   if (c.ns(ic) == 0)
      c.x(ic) = xsub(i) - dsubhx;
      c.y(ic) = ysub(i) + dsubhy;
   end;
   c.ns(ic) = c.ns(ic) + 1;
   c.ul(ic) = i;
   maxic = max([maxic, ic]);
   ic = findcorner(c, xsub(i)+dsubhx, ysub(i)+dsubhy, tol);
   if (c.ns(ic) == 0)
      c.x(ic) = xsub(i) + dsubhx;
      c.y(ic) = ysub(i) + dsubhy;
   end;
   c.ns(ic) = c.ns(ic) + 1;
   c.ur(ic) = i;
   maxic = max([maxic, ic]);
end;
c.x = c.x(1:maxic);
c.y = c.y(1:maxic);
c.ns = c.ns(1:maxic);
c.ll = c.ll(1:maxic);
c.lr = c.lr(1:maxic);
c.ul = c.ul(1:maxic);
c.ur = c.ur(1:maxic);
nc = maxic;
if (nargin > 2)
   if ((nc == prod(size(xact))) & (nc == prod(size(yact))))
      ind = zeros(nc,1);
      for ic=1:nc
         indi = find((abs(xact(ic)-c.x)+abs(yact(ic)-c.y)) < (2*tol));
         if (prod(size(indi)) == 1)
            ind(ic) = indi;
         else
            warning('Cannot find actuator.');
            return;
         end;
      end;
      c.x = c.x(ind);
      c.y = c.y(ind);
      c.ns = c.ns(ind);
      c.ll = c.ll(ind);
      c.lr = c.lr(ind);
      c.ul = c.ul(ind);
      c.ur = c.ur(ind);
   else
      warning('Number of actuators do not match');
      return;
   end;
end;
ptos = zeros(2*ns, nc);
for ic=1:nc
   if (c.ll(ic) > 0)
      ptos(     c.ll(ic),ic) = -0.5;
      ptos(ns + c.ll(ic),ic) = -0.5;
   end;
   if (c.lr(ic) > 0)
      ptos(     c.lr(ic),ic) =  0.5;
      ptos(ns + c.lr(ic),ic) = -0.5;
   end;
   if (c.ul(ic) > 0)
      ptos(     c.ul(ic),ic) = -0.5;
      ptos(ns + c.ul(ic),ic) =  0.5;
   end;
   if (c.ur(ic) > 0)
      ptos(     c.ur(ic),ic) =  0.5;
      ptos(ns + c.ur(ic),ic) =  0.5;
   end;
end;