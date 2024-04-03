rsx=[ones(nsub,1);zeros(nsub,1)];
rsy=[zeros(nsub,1);ones(nsub,1)];
if (size(mvtos,1) > (2*nsub))
   rsx = [rsx;zeros(nsub,1)];
   rsy = [rsy;zeros(nsub,1)];
end;
rb=pinv(mvtos);
rvx=rb*rsx;
rvy=rb*rsy;
rq=orth([rvx,rvy]);
rp=eye(size(rq,1))-(rq*rq');
rp1=eye(size(rsx,1))-((1/nsub)*(rsx*rsx' + rsy*rsy'));
recon=rp*rb*rp1';