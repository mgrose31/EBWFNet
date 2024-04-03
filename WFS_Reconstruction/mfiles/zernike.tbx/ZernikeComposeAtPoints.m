function [z]=ZernikeComposeAtPoints(v,x,y,zernikeRadius,varargin)
% SYNTAX:
% [z]=ZernikeComposeAtPoints(v,x,y,zernikeRadius,[aperture], ...
%           [normalization],[ordering])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab version of the WaveTrain ZernikeCompose system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% v [ ] = vector of Zernike coefficients (starting with tilts)
% x,y [] = vectors or arrays specifying locations of points for evaluation
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: ZernikeComposeAtPoints.m 3465 2012-05-16 17:24:24Z keith $

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

if numel(x) ~= numel(y)
    error('ZernikeComposeAtPoints: x and y must be same size!')
end
z = zeros(size(x));

for i=1:numel(x)
    r=sqrt(x(i)^2+y(i)^2)/zernikeRadius;
    theta=atan2(y(i),x(i));
    
    for ii=1:length(v);
        if(v(ii)~=0)
            [l,n]=TwoParameter(ii,ordering);
            z(i)=z(i)+v(ii)*RealZernike(l,n,r,theta,normalization,aperture);
        end
    end
    
end

return