[ru,rs,rv]=svd(mvtos);
rsv=diag(rs);
rmaxsv=max(rsv);
rnsv = size(rsv,1);
rsv1=ones(rnsv,1)./rsv;
if (exist('ns','var'))
   if (ns > 0), rsv1((rnsv-ns+1):rnsv)=zeros(ns,1); end;
   rsn1=[diag(rsv1),zeros(size(rsv1,1),size(ru,1)-size(rsv1,1))];
   recon=rv*rsn1*ru';
   fprintf(1,'The tilt-included reconstructor WAS computed.\n');
else
   fprintf(1,'The tilt-included reconstructor WAS NOT computed.\n');
end;
