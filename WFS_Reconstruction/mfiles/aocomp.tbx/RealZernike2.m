function [zz]=RealZernike(l,m,r,theta)
% SYNTAX: 
% [zz]=RealZernike(l,m,r,theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% l [ ] = 
% m [ ] = 
% r [ ] = 
% theta [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% zz [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: RealZernike2.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

luse=abs(l);
zz=sqrt((m+1)/pi).*RGZ(luse,m,r);
if (l<=0) 
    zz=zz.*cos(luse*theta);
else 
    zz=zz.*sin(luse*theta);
end
return;
