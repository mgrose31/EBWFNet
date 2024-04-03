function [z,xx,yy]=ZernikeCompose(v,nxy,dxy,zernikeRadius,varargin)
% SYNTAX:
% [z,xx,yy]=ZernikeCompose(v,nxy,dxy,zernikeRadius,[aperture], ...
%           [normalization],[ordering])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab version of the WaveTrain ZernikeCompose system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% v [ ] = vector of Zernike coefficients (starting with tilts)
% nxy,dxy [ ] = size and spacing of output grid (WaveTrain gwoom)
% zernikeRadius [ ] = radius of Zernike polynomials
% aperture(true) [ ] = only calculate within zernikeRadius (true/false)
% normalization(0) [ ] = normalization scheme
% ordering(0) [ ] = ordering scheme
%
%       Possible normalizations are:
%       0 - none
%       1 - overlap integral = 1/pi (Born & Wolf ?)
%       2 - phase variance = 1
%       3 - overlap integral = 1 (Noll, Malacara)
%
%       Possible orderings are:
%       0 - depreciated; do not use - you should get warnings
%       1 - Noll
%       2 - Malacara
%       3 - Wyant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% z [ ] = OPD due to Zernike Modes (same units as coefficients, v)
% xx [ ] = X grid coordinates matrix in form as output by meshgrid
% yy [ ] = Y grid coordinates matrix in form as output by meshgrid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: ZernikeCompose.m 3036 2010-09-23 21:22:50Z amoran $

%% BEGIN_CODE 

if (nargin>=5)
    aperture=varargin{1};
else
    aperture=true;
end;

if (nargin>=6)
    normalization=varargin{2};
else
    normalization=0;
end;

if (nargin>=7)
    ordering=varargin{3};
else
    ordering=0;
end;

if aperture
    % Clip grids to aperture size to minimize calculations below
    mxy = 2*ceil(zernikeRadius/dxy +1);
    mxy=min(mxy,nxy);% In case input grid was <= aperture
    x=(-mxy/2:1:mxy/2-1)*dxy/zernikeRadius; % Gwoom : Origin is at mxy/2 point with 0 based indexing
    z=zeros(mxy,mxy);
else
    x=(-nxy/2:1:nxy/2-1)*dxy/zernikeRadius; % Gwoom : Origin is at nxy/2 point with 0 based indexing
    z=zeros(nxy,nxy);
end

[xx,yy]=meshgrid(x);
r=sqrt(xx.^2+yy.^2);
theta=atan2(yy,xx);

for ii=1:length(v);
    if(v(ii)~=0)
        [l,n]=TwoParameter(ii,ordering);
        z=z+v(ii)*RealZernike(l,n,r,theta,normalization,aperture);
    end
end

if aperture
    % Pad with zeros to full size grid
    zp=zeros(nxy,nxy);
    zp((1:mxy)+(nxy-mxy)/2,(1:mxy)+(nxy-mxy)/2)=z; % Preserve origin at nx/2 (mx/2) with 1 based indexing
    z=zp;clear zp
    x=(-nxy/2:1:nxy/2-1)*dxy/zernikeRadius;
    [xx,yy]=meshgrid(x);
end

xx=xx*zernikeRadius;
yy=yy*zernikeRadius;

return