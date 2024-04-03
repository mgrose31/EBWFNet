function [apermask]=makemask(nacr,xcen,ycen,rin,rout);
% SYNTAX:
% [apermask]=makemask(nacr,xcen,ycen,rin,rout)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Generate (nacr x nacr) aperture mask for specified center, inner & outer 
% radius
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: makemask.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE
 
for ii=1:nacr;
  for jj=1:nacr;
    dist(ii,jj) = sqrt((xcen-ii)^2+(ycen-jj)^2);
  end
end
numaps=nacr^2;
reshape(dist,1,numaps);
idx=find((dist<rout)&(dist>rin));
apmask = zeros(1,numaps);
apmask(idx) = 1;
apermask = find(apmask>.9);
