function SRA = opd2slope(xsub, ysub, xopd, yopd)
% SYNTAX:
% SRA = opd2slope(xsub, ysub, xopd, yopd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Creates a DM displacment (OPD) to slope matrix Based on lsfptos; 
% modified from actuator grid to work with opd grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% xsub, ysub [ ] = The nsub vectors of subaperture locations.
% xopd, yopd [ ] = The nopd vectors of opd grid-point locations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% SRA [ ] = The 2 nsub x nopd matrix relating DM displacements to slopes.
% Notes:
%    The matrix will create slopes from mirror displacment.
%    Because of reflection, this is actually 1/2 of the OPD.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
% AUTHOR: Keith@MZA 3/30/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: opd2slope.m 3141 2010-11-11 21:06:28Z keith $

%% BEGIN_CODE

msg = nargchk(4,4,nargin);
if (~isempty(msg))
   warning(['opd2slope: ', msg]);
   help opd2slope;
   return;
end;
disp(['Calculating OPD to subaperture slope matrix ...'])
%
ns = length(xsub);
dsubs = abs(xsub(2:ns) - xsub(1:(ns-1)));
indnz = find(dsubs);
dsub = min(dsubs(indnz));
dsubhx = dsub / 2.0;
dsubs = abs(ysub(2:ns) - ysub(1:(ns-1)));
indnz = find(dsubs);
dsub = min(dsubs(indnz));
dsubhy = dsub / 2.0;
%
% Assume subaperture corners always coincide with some opd grid location
no = length(xopd);
dopds = abs(xopd(2:no) - xopd(1:(no-1)));
indnz = find(dopds);
dopd = min(dopds(indnz));
% Test that OPD spacing is an integer division of subap spacing
if abs(rem(dsub,dopd)/dopd) > 0.001 && abs(rem(dsub,dopd)/dopd) < 0.999
    error('opd2slope.m: OPD grid does not subsample subap grid!')
end
% Number of opd grid points to a side of the subaperture.
nos=round(dsub/dopd);
disp(['OPD grid is ',int2str(nos),' times finer than subaperture grid'])
%
% tol = 0.5 * 2*min([dsubhx,dsubhy])/nos; % Approx 1/2 OPD grid spacing
tol = 0.1 * 2*min([dsubhx,dsubhy])/nos; % Approx 1/10 OPD grid spacing
% List of grid points on corners & sides of subapertures
c.x = zeros(1,4*nos*ns);
c.y = c.x;
c.ns = c.x;
c.ll = c.x;
c.lr = c.x;
c.ul = c.x;
c.ur = c.x;
c.d = c.x; % down
c.u = c.x; % up
c.l = c.x; % left
c.r = c.x; % right
maxic = 0;
%disp(['Initial size of corner structure is ',int2str(size(c.x))])
%
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
   %
   % Find side points for this subaperture
   for j=1:nos-1
       % Lower edge
       xside=xsub(i)-dsubhx+j*dopd;
       yside=ysub(i)-dsubhy;
       ic = findcorner(c, xside, yside, tol);
       if (c.ns(ic) == 0)
           c.x(ic) = xside;
           c.y(ic) = yside;
       end;
       c.ns(ic) = c.ns(ic) + 1;
       c.d(ic) = i;
       maxic = max([maxic, ic]);
       % Upper edge
       xside=xsub(i)-dsubhx+j*dopd;
       yside=ysub(i)+dsubhy;
       ic = findcorner(c, xside, yside, tol);
       if (c.ns(ic) == 0)
           c.x(ic) = xside;
           c.y(ic) = yside;
       end;
       c.ns(ic) = c.ns(ic) + 1;
       c.u(ic) = i;
       maxic = max([maxic, ic]);
       % Left edge
       xside=xsub(i)-dsubhx;
       yside=ysub(i)-dsubhy+j*dopd;
       ic = findcorner(c, xside, yside, tol);
       if (c.ns(ic) == 0)
           c.x(ic) = xside;
           c.y(ic) = yside;
       end;
       c.ns(ic) = c.ns(ic) + 1;
       c.l(ic) = i;
       maxic = max([maxic, ic]);
       % Right edge
       xside=xsub(i)+dsubhx;
       yside=ysub(i)-dsubhy+j*dopd;
       ic = findcorner(c, xside, yside, tol);
       if (c.ns(ic) == 0)
           c.x(ic) = xside;
           c.y(ic) = yside;
       end;
       c.ns(ic) = c.ns(ic) + 1;
       c.r(ic) = i;
       maxic = max([maxic, ic]);
   end
end;
c.x = c.x(1:maxic);
c.y = c.y(1:maxic);
c.ns = c.ns(1:maxic);
c.ll = c.ll(1:maxic);
c.lr = c.lr(1:maxic);
c.ul = c.ul(1:maxic);
c.ur = c.ur(1:maxic);
c.d = c.d(1:maxic);
c.u = c.u(1:maxic);
c.l = c.l(1:maxic);
c.r = c.r(1:maxic);
nc = maxic;
%disp(['Final size of corner structure is ',int2str(size(c.x))])
%
% Find & index OPD points at corners & sides
for ic=1:nc
    indi = find((abs(xopd-c.x(ic)) < tol) & (abs(yopd-c.y(ic)) < tol));
    if (prod(size(indi)) == 1)
        c.ind(ic) = indi;
    elseif(prod(size(indi)) == 0)
        warning(['Cannot find opd point at (',num2str(c.x(ic)),',',num2str(c.y(ic)),')']);
        return;
    else
        warning(['Multiple opd points near (',num2str(c.x(ic)),',',num2str(c.y(ic)),')']);
        disp('   X       Y')
        disp([xopd(indi) yopd(indi)])
        return;
    end;
end;
%
SRA = zeros(2*ns, no);
for ic=1:nc
   if (c.ll(ic) > 0)
      SRA(     c.ll(ic),c.ind(ic)) = -0.5;
      SRA(ns + c.ll(ic),c.ind(ic)) = -0.5;
   end;
   if (c.lr(ic) > 0)
      SRA(     c.lr(ic),c.ind(ic)) =  0.5;
      SRA(ns + c.lr(ic),c.ind(ic)) = -0.5;
   end;
   if (c.ul(ic) > 0)
      SRA(     c.ul(ic),c.ind(ic)) = -0.5;
      SRA(ns + c.ul(ic),c.ind(ic)) =  0.5;
   end;
   if (c.ur(ic) > 0)
      SRA(     c.ur(ic),c.ind(ic)) =  0.5;
      SRA(ns + c.ur(ic),c.ind(ic)) =  0.5;
   end;
   if (c.d(ic) > 0)
      SRA(ns + c.d(ic),c.ind(ic)) = -1.0;
   end;
   if (c.u(ic) > 0)
      SRA(ns + c.u(ic),c.ind(ic)) =  1.0;
   end;
   if (c.l(ic) > 0)
      SRA(     c.l(ic),c.ind(ic)) = -1.0;
   end;
   if (c.r(ic) > 0)
      SRA(     c.r(ic),c.ind(ic)) =  1.0;
   end;
end;
%
% Factor of two due to reflection
% 1/dsub for subaperture width
% 1/nos to scale/average contributions from all opd points
SRA = 2*SRA/(dsub*nos);
disp(['Calculated OPD to subaperture slope matrix'])
