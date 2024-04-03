[ru,rs,rv]=svd(mvtos);
rsx=[ones(nsub,1);zeros(nsub,1)];
rsy=[zeros(nsub,1);ones(nsub,1)];
if (size(ru,1) > (2*nsub))
   rsx = [rsx;zeros(nsub,1)];
   rsy = [rsy;zeros(nsub,1)];
end;
rb=pinv(mvtos);
imasact=find(actype>0);
rvx=yact(imasact);%rvx=rb*rsx;
rvy=xact(imasact);%rvy=rb*rsy;
rq=orth([rvx,rvy]);
rp=eye(size(rq,1))-(rq*rq');
rp1=eye(size(rsx,1))-((1/nsub)*(rsx*rsx' + rsy*rsy'));
rsv=diag(rs);
rmaxsv=max(rsv);
rnsv = size(rsv,1);
rsv1=ones(rnsv,1)./rsv;
if (exist('ns','var'))
   if (ns > 0), rsv1((rnsv-ns+1):rnsv)=zeros(ns,1); end;
   rsn1=[diag(rsv1),zeros(size(rsv1,1),size(ru,1)-size(rsv1,1))];
   rb1=rv*rsn1*ru';
   recon=rp*rb1*rp1';
   fprintf(1,'The reconstructor WAS computed.\n');
else
   fprintf(1,'The reconstructor WAS NOT computed.\n');
end;
