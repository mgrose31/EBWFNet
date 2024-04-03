function [zz] = SigmaGZ(l,m,r)
% SYNTAX: 
% [zz] = SigmaGZ(l,m,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: SigmaGZ.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

l = abs(l);
zz = r.*0.0;
sstop = (m-l)/2.0;
if l== 0
    return;
end
for s=0:sstop;
    coef=1.0/(((m+l)/2-s+1)*((m-l)/2-s+1));
    power1=l-1;
    power2=m-2*s+1;
    acoef=aGZ(l,m,s);
    totalcoef=acoef*coef/(4.0*l);
    zz=zz+totalcoef*((m-2*s+2)*r.^power1-l*r.^power2);
end
return;