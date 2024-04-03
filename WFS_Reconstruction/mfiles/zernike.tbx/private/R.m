function [zz,rcoefs]=R(n,l,r)
% SYNTAX:
% [zz,rcoefs]=R(n,l,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: R.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE 

m=(n-l)/2;
sstop=m;
zz=(0.0).*r;
for s=0:sstop
    pow=n-2*s;
    coef=a(l,n,s);
    switch pow
        case 0
            zz=zz+(coef.*ones(size(r)));
        case 1
            zz=zz+(coef.*r);
        case 2
            zz=zz+(coef.*(r.*r));
        otherwise
            zz=zz+(coef.*power(r,pow));
    end
    rcoefs(pow+1) = coef;
end;
return;