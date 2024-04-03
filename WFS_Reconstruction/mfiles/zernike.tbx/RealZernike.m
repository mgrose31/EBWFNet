function [zz]=RealZernike(l,n,r,theta,varargin)
% SYNTAX:
% [zz]=RealZernike(l,n,r,theta,[normalization],[aperture])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Matlab version of WaveTrain rZernike() in processing.lib/ZernikeRoutines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% n [ ] = radial order (starting at 1 for tilts)
% l [ ] = angular order (one of -n:2:n)
% r [ ] = spherical coordinates: radius
% theta [ ] = spherical coordinates: angle
% normalization(0) [ ] = normalization scheme
% aperture(true) [ ] = only calculate within zernikeRadius (true/false)
%
%       Possible normalizations are:
%       0 - none
%       1 - overlap integral = 1/pi (Born & Wolf ?)
%       2 - phase variance = 1
%       3 - overlap integral = 1 (Noll, Malacara)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% zz [ ] = OPD due to Zernike Modes (dimensionless)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: RealZernike.m 3141 2010-11-11 21:06:28Z keith $

%% BEGIN_CODE 

if (nargin>=5)
    normalization=varargin{1};
else
    normalization=0;
end;
if (nargin>=6)
    aperture=varargin{2};
else
    aperture=true;
end;

if (n==0)
    zz=ones(size(r));
else
    zz= R(n,abs(l),r);

    if (l>0)
        zz=zz.*sin(l*theta);
    else
        zz=zz.*cos(l*theta);
    end
end

if (normalization)
    zz=zz*ZernikeNormalizationFactor(l,n,normalization);
end;

if (aperture); zz=zz.*(r<=1.0); end;

return

