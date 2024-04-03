function [c]=mknollv(nzer,rad);
% SYNTAX:
% [c]=mknollv(nzer,rad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% nzer [ ] = 
% rad [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% c [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: mknollv.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE
 
l=5/3;
fac1 = -3.44*pi*sqrt(pi)*gamma(11/3)/gamma(-5./6.)/gamma(7/3)*gamma(14/3)/(2^(13/3))*rad^2;
[rws,cls]=meshgrid([1:1:nzer]);
c = zeros(1,prod(size(rws)));
sig = c;
m = c;
rws = reshape(rws,1,prod(size(rws)));
cls = reshape(cls,1,prod(size(cls)));
[np,mi,sincosi]=pzrntonmv(rws);
[nq,mj,sincosj]=pzrntonmv(cls);
yy = find((mj == mi) & (sincosi == sincosj));
m(yy) = mj(yy);
sig = (-1*ones(1,prod(size(rws)))).^((np-m+nq-m)/2);
fac = sig*fac1;
c(yy) = fac(yy).*sqrt((np(yy) + 1) .* (nq(yy) + 1)).*gamma((np(yy) + nq(yy))/2.0 - 5.0 / 6.0)...
        ./gamma((np(yy) - nq(yy)) / 2.0 + 17.0 / 6.0)./gamma((nq(yy)- np(yy)) / 2.0 + 17.0 / 6.0)...
        ./gamma((np(yy) + nq(yy)) / 2.0 + 23.0 / 6.0);
c=reshape(c,nzer,nzer);
