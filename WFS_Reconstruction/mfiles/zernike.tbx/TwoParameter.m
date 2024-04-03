function [l,n]=TwoParameter(index,ordering)
% SYNTAX:
% [l,n]=TwoParameter(index,ordering)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Matlab version of zern_ToNandL() in WaveTrain processing.lib/ZernikeRoutines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% index [ ] = 1D Zernike index (starting with tilts)
% ordering [ ] = ordering scheme
%
%       Possible orderings are:
%       0 - depreciated; do not use - you should get warnings
%       1 - Noll
%       2 - Malacara
%       3 - Wyant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUPTUTS:
% n [ ] = radial order (starting at 1 for tilts)
% l [ ] = angular order (one of -n:2:n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: TwoParameter.m 3141 2010-11-11 21:06:28Z keith $

%% BEGIN_CODE 

if ordering == 1
    [l,n]=Noll(index); % Noll ordering
elseif ordering == 2
    [l,n]=Malacara(index); % Malacara ordering
elseif ordering == 3
    [l,n]=Wyant(index); % Wyant ordering
else
    warning('MATLAB:ZernikeOrdering','Using deprecated ordering scheme; use a non-zero value for the *ordering* parameter')
    [l,n]=Old_Noll(index);
end

return;