function [val]=pzernikev(n,m,r);
% SYNTAX:
% [val]=pzernikev(n,m,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: pzernikev.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE
 
val=ones(1,length(r));
z = val;

if m == 0
   tmp=find(r==0);
   z(tmp) = 1;
else
   z = r.^m;
end

if n==m
   val = z;
else
   zm = z;
   rsq = r.*r;
   mp2 = m+2;
   z = (mp2*rsq-m-1).*zm;
   if n==mp2
      val=z;
   else
      for n1=mp2:2:n-2
          n2 = n1+2;
          azp=n1*(n2-m)*(n2+m);
          azm=-1*n2*(n1+m)*(n1-m);
          az = 4.*n1*n2*(n1+1)*rsq-n1*(n2-m)^2-n2*(n1+m)^2;
          zp = (az.*z+azm*zm)/azp;
          zm = z;
          z = zp;
          val = z;
      end
   end
end
