%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: TestZernikeFit.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE 

close all;clear all;

nxy = 512;
dxy = 7e-6;
r = 8.0e-4;
xy = dxy*([0 nxy-1]-nxy/2);
nzern = 24;           % Number of Zernike coefficients to fit

apFctn = aperture(nxy,dxy,r);

useAperture = true;
normalization = 2;% 0 - none, 1 -overlap=1/pi, 2 - pv=1, 3 - Noll/Malaraca
ordering=2;       % 1 -Noll, 2 - Malacara

for ii = 1:20
    v=zeros(1,ii);v(ii)=1;
    phase = ZernikeCompose(v,nxy,dxy,r,useAperture,normalization,ordering); % ith Zernike mode

    v=ZernikeDecompose(phase,dxy,dxy,nzern,r,normalization,ordering);   % Zernike decomposition of phase
    phaseFit=ZernikeCompose(v,nxy,dxy,r,useAperture,normalization,ordering); % Retrieved phase
    [l,n]=TwoParameter(ii,ordering);
    figure('name',ZernikeName(l,n));
    clim=[min(phase(:)) max(phase(:))];
    subplot(2,2,1);imagesc(xy,xy,phase,clim);colorbar;axis image;axis([-1 1 -1 1]*r);title('Phase')
    subplot(2,2,2);imagesc(xy,xy,phaseFit,clim);colorbar;axis image;axis([-1 1 -1 1]*r);title('Phase fit')
    subplot(2,2,3);imagesc(xy,xy,phase-phaseFit,clim);colorbar;axis image;axis([-1 1 -1 1]*r);title('Phase error')
    subplot(2,2,4);imagesc(xy,xy,phase-phaseFit);colorbar;axis image;axis([-1 1 -1 1]*r);title('Phase error')
end