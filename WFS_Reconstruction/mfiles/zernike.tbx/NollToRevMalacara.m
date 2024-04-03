function [vM] = NollToRevMalacara(vN)
% SYNTAX:
% [vM] = NollToRevMalacara(vN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS: 
% vN [ ] = Zernike coefficients in NOLL order, starting with tilts (no piston)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS: 
% vM [ ] = Zernike coefficients in reverse Malacara order, starting with 
%          tilts (no piston)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: NollToRevMalacara.m 3036 2010-09-23 21:22:50Z amoran $

%% BEGIN_CODE 

vM = zeros(1,2*length(vN));
jMax = 0;
for i=1:length(vN)
    [l,n]=Noll(i);
    %m=(n-l)/2;% Malacara (13.21)
    m=(n+l)/2;% Malacara (13.21)
    j=n*(n+1)/2 + m + 1; % Malacara (13.25)
    j =j-1;
    vM(j) = sqrt(pi)*vN(i);
    if j>jMax
        jMax=j;
    end
end
vM=vM(1:jMax);
