function [vN] = RevMalacaraToNoll(vM)
% SYNTAX:
% [vN] = RevMalacaraToNoll(vM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% vM [ ] = Zernike coefficients in reverse Malacara order, starting with 
%          tilts (no piston)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% vN [ ] = Zernike coefficients in NOLL order, starting with tilts 
%          (no piston)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: RevMalacaraToNoll.m 3036 2010-09-23 21:22:50Z amoran $

%% BEGIN_CODE 

vN = zeros(1,2*length(vM));
jMax = 0;
for i=1:length(vM)
    [l,n]=Malacara(i);
    l=-l; % Reverse azimuthal order
    if l==0
        j=n*(n+1)/2+1;
    else
        j=n*(n+1)/2+abs(l);
        % l -ve => cos() => j even
        if l < 0 && mod(j,2)
            j=j+1;
        elseif l > 0 && ~mod(j,2)
            j=j+1;
        end
    end
    j=j-1;
    vN(j) = vM(i)/sqrt(pi);
    if j>jMax
        jMax=j;
    end
end
vN=vN(1:jMax);
