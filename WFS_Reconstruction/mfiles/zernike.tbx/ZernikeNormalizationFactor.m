function [x]=ZernikeNormalizationFactor(l,n,normalization)
% SYNTAX: 
% [x]=ZernikeNormalizationFactor(l,n,normalization)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Matlab version of WaveTrain zernNorm() in WaveTrain 
% processing.lib/ZernikeRoutines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% n [ ] = radial order (starting at 1 for tilts)
% l [ ] = angular order (one of -n:2:n)
% normalization [ ] = normalization scheme
%
%       Possible normalizations are:
%       0 - none
%       1 - overlap integral = 1/pi (Born & Wolf ?)
%       2 - phase variance = 1
%       3 - overlap integral = 1 (Noll, Malacara)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% x [ ] = normalization factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: ZernikeNormalizationFactor.m 3141 2010-11-11 21:06:28Z keith $

%% BEGIN_CODE 

x=1; % No normalization by default

if normalization == 1 % Overlap integral = 1/pi (Born & Wolf ?)
    x=sqrt(2.0*(n+1.0)/pi);
    if (mod(n,2)==0 && l==0); x=x./sqrt(2); end;
elseif normalization == 2 % Phase Variance = 1
    x=0.5;
    if (mod(n,4)==0 && l==0) x=1/sqrt(2); end;
	if (n==4 && l==0) x=1/1.5; end;
    if (n==8 && l==0) x=1/(1+3/7); end;
elseif normalization == 3 % Overlap integral = 1 (Noll/Malacara)
        x=sqrt(2.0*(n+1.0));
    if (mod(n,2)==0 && l==0); x=x./sqrt(2); end;
end

return;