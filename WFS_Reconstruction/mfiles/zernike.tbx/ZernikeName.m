function [name]=ZernikeName(l,n)
% SYNTAX:
% [name]=ZernikeName(l,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Returns the common name of the Zernke term using 2-parameter indexing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% n [ ] = radial order (starting at 1 for tilts)
% l [ ] = angular order (one of -n:2:n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% name [ ] = common name of Zernike mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: ZernikeName.m 3141 2010-11-11 21:06:28Z keith $

%% BEGIN_CODE 

if mod(n-l,2) || abs(l)>n
    name='Bad (l,n) parameters';
elseif l==-1 && n==1
    name='X Tilt';
elseif l==1 && n==1
    name='Y Tilt';
elseif l==0 && n==2
    name='Focus';
elseif l==2 && n==2
    name='+/-45 Astigmatism';
elseif l==-2 && n==2
    name='0/90 Astigmatism';
elseif l==1 && n==3
    name='Y Coma';
elseif l==-1 && n==3
    name='X Coma';
elseif l==3 && n==3
    name='Y Trefoil';
elseif l==-3 && n==3
    name='X Trefoil';
elseif l==0 && n==4
    name='Sph. Aberr.';
elseif l==-2 && n==4
    name='0/90 Secondary astigmatism';
elseif l==2 && n==4
    name='+/-45 Secondary astigmatism';
elseif l==-4 && n==4
    name='0/90 Quadrafoil';
elseif l==4 && n==4
    name='+/-45 Quadrafoil';
elseif l==-1 && n==5
    name='X Secondary coma';
elseif l==1 && n==5
    name='Y Secondary coma';
elseif l==-3 && n==5
    name='X Secondary trefoil';
elseif l==3 && n==5
    name='Y Secondary trefoil';
elseif l==-5 && n==5
    name='X Pentafoil';
elseif l==5 && n==5
    name='Y Pentafoil';
else
%    name=sprintf('Z(%i,%i)',n,l);
    name=sprintf('Z_{%i}^{%i}',n,l);
end

return;