function [zz]=dZdr(l,m,r)
zz=(0.0).*r;
l=abs(l);
%
sstop=(m-l)/2;
for s=0:sstop
    power=m-2*s;
    zz=zz+(aGZ(l,m,s)*power.*r.^(power-1));
end;
%
return;