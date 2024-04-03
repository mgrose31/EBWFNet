function [v]=ZernikeDecompose(z,dx,dy,zernikeCount,zernikeRadius,varargin)
% SYNTAX:
% [v]=ZernikeDecompose(z,dx,dy,zernikeCount,zernikeRadius, ...
%     [normalization],[ordering],[ls_tol])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Matlab version of the WaveTrain ZernikeDecompose system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% z [ ] = array of phase or OPD (WaveTrain gwoom)
% dx,dy [ ] = spacing of input grid 
% zernikeCount [ ] = total number of Zernike coefficients to extract (starting with tilts)
% zernikeRadius [ ] = radius of Zernike polynomials
% normalization(0) [ ] = normalization scheme
% ordering(0) [ ] = ordering scheme
% ls_tol(1) [] = tolerance of non-orthogonality before switching to least squares
%                   1=never do least squares, 0=always do least squares,
%                   0<ls_tol<1 = switch to least squares when that fraction
%                   of pixels are obscured (by annulus, etc)
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
% v [ ] = vector of Zernike coefficients (starting with tilts, same units 
%         as OPD, z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: ZernikeDecompose.m 3141 2010-11-11 21:06:28Z keith $

%% BEGIN_CODE 

if (nargin>=6)
    normalization=varargin{1};
else
    normalization=0;
end;

if (nargin>=7)
    ordering=varargin{2};
else
    ordering=0;
end;

if (nargin>=8)
    ls_tol=varargin{3};
else
    ls_tol=1;
end;

[ny,nx]=size(z);

% x=(-nx/2:1:nx/2-1)*dx/zernikeRadius; % Gwoom : Origin is at nx/2 point with 0 based indexing
% y=(-ny/2:1:ny/2-1)*dy/zernikeRadius;

% Clip grids to aperture size to minimize calculations below
mx = 2*ceil(zernikeRadius/dx +1);
my = 2*ceil(zernikeRadius/dy +1);
if mx > nx || my > ny
    warning('ZernikeDecompose: input mesh is smaller than Zernike diameter');
    mx=min(nx,mx);
    my=min(ny,my);
    orthogonal = false;% Can't be orthogonal
else
    orthogonal = true; % At least for now
end
x=(-mx/2:1:mx/2-1)*dx/zernikeRadius; % Gwoom : Origin is at mx/2 point with 0 based indexing
y=(-my/2:1:my/2-1)*dy/zernikeRadius;
z=z((1:my)+(ny-my)/2,(1:mx)+(nx-mx)/2); % Preserve origin at nx/2 (mx/2) with 1 based indexing

[xx,yy]=meshgrid(x,y);
r=sqrt(xx.^2 + yy.^2);
theta=atan2(yy,xx);

% Determine if Zernikes are at least approximately orthogonal
npoints =sum(r(:)<=1); % Number of points within radius
nNaNs = sum(~isfinite(z(r<=1))); % Number of unknown points within radius
ratNaN = nNaNs/npoints;
orthogonal = orthogonal && ratNaN<ls_tol; % Test for approximate orthogonality
clear npoints nNaNs ratNaN
%
if ls_tol >=1
    if ~orthogonal
        warning('Zernike Decompose using integral calculation but Zernikes do not appear to be orthogonal!')
    end
    orthogonal = true;
elseif ls_tol <=0
    orthogonal = false;
else
    if ~orthogonal
        disp('Zernike Decompose switching to least squares fit (non-orthogonal)')
    end
end

% Normalize sum by number of points in circle
% Also ignore NaN values in phase sheet
mask = (r<=1).*isfinite(z);
npoints =sum(mask(:));
z(~isfinite(z))=0;

if orthogonal
    % Orthonormal; integral fit
    z = z - mean(z(mask>0)); % Remove any large piston term that could screw up not-quite-orthogonal fits
    v=zeros(1,zernikeCount);
    for ii=1:zernikeCount
        [l,n]=TwoParameter(ii,ordering);
        zern=RealZernike(l,n,r,theta,3,true); % Flags normalize integral to 1 on unit circle and apply aperture
        v(ii)=sum(zern(:).*z(:))/npoints; % Integral over unit circle
        % Scale coefficients to required normalization
        v(ii) = v(ii) * ZernikeNormalizationFactor(l,n,3)/ZernikeNormalizationFactor(l,n,normalization);
    end
else
    % Non-orthogonal; least-squares fit
    r = r(mask>0);
    theta = theta(mask>0);
    z = z(mask>0);
    X=NaN*zeros(numel(z),1+zernikeCount);
    X(:,1)=1; % Piston
    for ii=1:zernikeCount
        [l,n]=TwoParameter(ii,ordering);
        zern=RealZernike(l,n,r,theta,3,true); % Flags normalize integral to 1 on unit circle and apply aperture
        X(:,ii+1)=zern(:) / (ZernikeNormalizationFactor(l,n,3)/ZernikeNormalizationFactor(l,n,normalization));
    end
    v=(X\z(:))';
    v=v(2:end); % Remove piston
end

return