function [zz]=RGZ(l,m,r)
% SYNTAX:
% [zz]=RGZ(l,m,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: RGZ.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

l = abs(l);
sstop=(m-l)/2;
zz=(0.0).*r;
for s=0:sstop
    power=m-2*s;
    zz=zz+(aGZ(l,m,s).*r.^(power));
end;
return;