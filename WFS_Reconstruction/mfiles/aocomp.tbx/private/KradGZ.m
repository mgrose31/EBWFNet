function [zz]=KradGZ(l,m,r)
% SYNTAX:
% [zz]=KradGZ(l,m,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: KradGZ.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

luse = abs(l);
zz=r.*0.0;
if (l==0)
    sstop=m/2;
    for s=0:sstop
        coef=aGZ(0,m,s)/(m-2*s+2);
        power1=m-2*s+1;
        zz=zz-coef*r.^power1;
    end
else   
    sstop=(m-luse)/2.0;
	for s=0:sstop;
        coef=(m-2*s+2)/(((m+l)/2-s+1)*((m-l)/2-s+1));
        power1=l-1;
        power2=m-2*s+1;
        acoef=aGZ(luse,m,s);
        totalcoef=acoef*coef/4.0;
        zz=zz+totalcoef*(r.^power1-r.^power2);
	end;
end
return;