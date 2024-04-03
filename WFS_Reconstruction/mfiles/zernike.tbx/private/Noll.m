function [l,n]=Noll(index)
% SYNTAX:
% [l,n]=Noll(index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: Noll.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE 

%function [l,n]=Noll(index)
% R. J. Noll, "Zernike polynomials and atmospheric turbulence," J. Opt. Soc. Am. 66, 207- (1976)
% Inputs:
% index - 1D Zernike index (starting with tilts)
%
% Outputs:
% n - radial order (starting at 1 for tilts)
% l - angular order (one of -n:2:n)

r=index+1; % Start with tilt, not piston
n=ceil((-3 + sqrt(9 + 8*(r-1)))/2); % From Malacara
if ~mod(n,2) % l for even n starts at 0
    absl=r-n*(n+1)/2-1;
    if mod(absl,2) % Add 1 to ever other value for steps of 2
        absl=absl+1;
    end
else % l for odd n starts at 1
    absl=r-n*(n+1)/2;
    if ~mod(absl,2)
        absl=absl-1; % Subtract 1 from ever other value for steps of 2
    end
end
if ~mod(r,2)
    l=-absl; % Even terms are cos, -ve l
else
    l=absl; % Odd terms are sin, +ve l
end

return;