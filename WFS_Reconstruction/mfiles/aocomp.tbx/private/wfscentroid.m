function [centroid, peakpix, maxval, pixsum] = wfscentroid(frame, wg, ...
                                            relthresh, thresh, type, mask)
% SYNTAX: 
% [centroid, peakpix, maxval, pixsum] = wfscentroid(frame, wg, relthresh, 
%                                       thresh, type, mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: wfscentroid.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

if (~exist('relthresh','var')), relthresh = 0; end;
if (isempty(relthresh)), relthresh = 0; end;
if (~exist('thresh','var')), thresh = 0; end;
if (isempty(thresh)), thresh = 0; end;
if (~exist('type','var')), type = 1; end;
if (isempty(type)), type = 1; end;
nx = wg.nx;
ny = wg.ny;
if (~exist('mask','var')), mask = ones(wg.nx, wg.ny); end;
if (isempty(mask)), mask = ones(wg.nx, wg.ny); end;
sm = size(mask);
if ((sm(1) == wg.xmax) & (sm(2) == wg.ymax))
   mask = slmtosamask(wg, mask);
end;
centroid = zeros(2,nx,ny);
peakpix = zeros(2,nx,ny);
maxval = zeros(nx,ny);
pixsum = zeros(nx,ny);
for jj=1:nx
   for ii=1:ny
      if (mask(jj,ii))
         [c, pp, mv] = sacentroid(frame(wg.u(jj):wg.d(jj),wg.l(ii):wg.r(ii)), relthresh, thresh, type);
         centroid(1,ii,jj) = (wg.l(ii) + c(2) - 1) - wg.y(ii); 
         centroid(2,ii,jj) = (wg.u(jj) + c(1) - 1) - wg.x(jj);
         peakpix(1,ii,jj) = (wg.l(ii) + pp(2) - 1) - wg.y(ii); 
         peakpix(2,ii,jj) = (wg.u(jj) + pp(1) - 1) - wg.x(jj);
         maxval(ii,jj) = mv;
         pixsum(ii,jj) = sum(sum(frame(wg.l(ii):wg.r(ii),wg.u(jj):wg.d(jj))));
      end;
   end;
end;