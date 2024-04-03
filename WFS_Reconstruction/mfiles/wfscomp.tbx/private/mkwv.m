function [wmat]=mkwv(procInfo,nz,nsimp);
% SYNTAX:
% [wmat]=mkwv(procInfo,nz,nsimp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computes the zernike coefficient to slopes matrix w. procInfo is a wfs 
% data array of structures nz the number of zernike coeeficients to be 
% considered and nsimp the number of steps in the simpson integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: mkwv.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE
 
xc = procInfo.xcenter;
yc = procInfo.ycenter;
radius= procInfo.outradius;
h=procInfo.spacing;
xi=procInfo.xsub;
yi=procInfo.ysub;
xi=xi.';
yi=yi.';
rsq=radius^2;
nsub=length(xi);
w=zeros(2*nsub,nz);
ind = [1:1:nsub];

hp = h/radius;
hpsq = hp*hp;
for jj=1:1:nsub
    xip(jj) = (xi(jj)-xc)/radius;
    yip(jj) = (yi(jj)-yc)/radius;
end

delh = hp/nsimp;
delhd3 = delh/3.0;
wsimp(1) = delhd3;
for jj=2:2:nsimp
        wsimp(jj) = 4.0*delhd3;
        wsimp(jj+1) = 2.0*delhd3;
end
wsimp(nsimp+1) = delhd3;

[n,m,sincos] = pzrntonmv([1:1:nz]);
cnp = zeros(1,length(n));
a = find(m==0);
cnp = (sqrt(2*(n+1)))/pi;
cnp(a) = (sqrt(1+n(a)))/pi;

xleft = ones(nsimp+1,1)*xip;
ylow = ones(nsimp+1,1)*yip;
xright = xleft+hp;
ytop = ylow+hp;

simpson = [0:1:nsimp]*delh;
simpsonb = simpson.'*ones(1,nsub);

x1b = xleft + simpsonb; 
x2b = ylow + simpsonb;

theta1b = atan2(ylow,x1b);
rho1b = sqrt(x1b.*x1b+ylow.*ylow);
theta2b = atan2(x2b,xright);
rho2b = sqrt(x2b.*x2b + xright.*xright);
theta3b = atan2(ytop,x1b);
rho3b = sqrt(ytop.*ytop+x1b.*x1b);       
theta4b = atan2(x2b,xleft);
rho4b = sqrt(x2b.*x2b + xleft.*xleft);

for k=1:1:nz
    if sincos(k) ==1
       tfun1 = sin(m(k)*theta1b)*cnp(k)/hpsq;
       tfun2 = sin(m(k)*theta2b)*cnp(k)/hpsq;
       tfun3 = sin(m(k)*theta3b)*cnp(k)/hpsq;
       tfun4 = sin(m(k)*theta4b)*cnp(k)/hpsq;
    else
       tfun1 = cos(m(k)*theta1b)*cnp(k)/hpsq;
       tfun2 = cos(m(k)*theta2b)*cnp(k)/hpsq;
       tfun3 = cos(m(k)*theta3b)*cnp(k)/hpsq;
       tfun4 = cos(m(k)*theta4b)*cnp(k)/hpsq;
    end

    zfun1 = pzernikev(n(k),m(k),reshape(rho1b,1,prod(size(rho1b))));
    z1 = reshape(zfun1,nsimp+1,nsub).*tfun1;
    z1v = wsimp*z1;
    z1s = z1v;

    zfun2 = pzernikev(n(k),m(k),reshape(rho2b,1,prod(size(rho2b))));
    z2 = reshape(zfun2,nsimp+1,nsub).*tfun2;
    z2v = wsimp*z2;
    z2s = z2v;

    zfun3 = pzernikev(n(k),m(k),reshape(rho3b,1,prod(size(rho3b))));
    z3 = reshape(zfun3,nsimp+1,nsub).*tfun3;
    z3v = wsimp*z3;
    z3s = z3v;

    zfun4 = pzernikev(n(k),m(k),reshape(rho4b,1,prod(size(rho4b))));
    z4 = reshape(zfun4,nsimp+1,nsub).*tfun4;
    z4v = wsimp*z4;
    z4s = z4v;

    w(ind,k) = (z2s-z4s).';
    w(ind+nsub,k) = (z3s-z1s).';
end
wmat = 1/rsq*w;
