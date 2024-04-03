function [dxxa,dyya,dxya]=cpaperstrt(nsub);
% SYNTAX:
% [dxxa,dyya,dxya]=cpaperstrt(nsub)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Produces the values of the x_self slope structure function for  nsub 
% displacements of the subapertures in the x-direction dxxa, in the 
% y-direction dyya, and along the diagonal direction dxya.  Recall that 
% the x-self slope structure function is a function of separation -- to 
% generate the graph of the function multiple the aperture separation by 
% the width of the sub-apertures for the graphs along the x & y directions 
% and by square root of two times the subaperture widths for the 
% diplacement along the diagonal.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: cpaperstrt.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE
 
d=1.;
global D;
D=d;

nsteps = 100;  % the number of intergation steps required
delx=1/nsteps;
xpos = [1:1:nsteps]*delx;
tmp = ((sqrt(1+xpos.^2)).^(5/3)-(abs(xpos)).^(5/3)).*tri(xpos);
term1 = 2*2*sum(tmp)*delx;  % extra factor of 2 takes into account the symmetry of the integrand 
factor = 6.88*(d^(-1/3));

for nn=1:1:nsub
    xn = 0;
    xk = (nn-1)*d;
    yn = 0;
    yk =0;
    xright = [1:1:nsteps]*delx + (yn-yk)/d;
    xleft = (yn-yk)/d - [1:1:nsteps]*delx;
    xsq = xright.*xright;
    t1 = (xn-xk)/d;
    pp1 = (sqrt((t1-1)^2+xsq)).^(5/3);
    pp2 = (sqrt((t1+1)^2+xsq)).^(5/3);
    pp3 = (sqrt((t1)^2+xsq)).^(5/3);
    tmp = (pp1+pp2-2*pp3).*tri(xright-((yn-yk)/d));
    tmp1 = (pp1+pp2-2*pp3).*tri(xleft-((yn-yk)/d));
    dxxa(nn)=(sum(tmp)+sum(tmp1))*delx;
end
dxxa=factor*(-1*dxxa+term1);
%
% compute the x-self slope structure function in the y-direction
%
for nn=1:1:nsub
   yn = -15*d; 
   yk = -15*d+(nn-1)*d;
   xn = 0;
   xk =0;
   xright = [1:1:nsteps]*delx + (yn-yk)/d;
   xleft = (yn-yk)/d - [1:1:nsteps]*delx;
   xsq = xright.*xright;
   t1 = (xn-xk)/d;
   pp1 = (sqrt((t1-1)^2+xsq)).^(5/3);
   pp2 = (sqrt((t1+1)^2+xsq)).^(5/3);
   pp3 = (sqrt((t1)^2+xsq)).^(5/3);
   tmp = (pp1+pp2-2*pp3).*tri(xright-((yn-yk)/d));
   tmp1 = (pp1+pp2-2*pp3).*tri(xleft-((yn-yk)/d));
   dyya(nn)=(sum(tmp)+sum(tmp1))*delx;
end
dyya=factor*(-1*dyya+term1);
%
% compute the x-self slope structure function along the diagonal
%
for nn=1:1:nsub
   xn = 0;
   xk = (nn-1)*d;
   yn = 0;
   yk = (nn-1)*d;
   xright = [1:1:nsteps]*delx + (yn-yk)/d;
   xleft = (yn-yk)/d - [1:1:nsteps]*delx;
   xsq = xright.*xright;
   t1 = (xn-xk)/d;
   pp1 = (sqrt((t1-1)^2+xsq)).^(5/3);
   pp2 = (sqrt((t1+1)^2+xsq)).^(5/3);
   pp3 = (sqrt((t1)^2+xsq)).^(5/3);
   tmp = (pp1+pp2-2*pp3).*tri(xright-((yn-yk)/d));
   tmp1 = (pp1+pp2-2*pp3).*tri(xleft-((yn-yk)/d));
   dxya(nn)=(sum(tmp)+sum(tmp1))*delx;
end
dxya=factor*(-1*dxya+term1);
%
%
function [y]=tri(x);
global D;
d=D;
y=0*x;
aa=find(abs(x) < 1.0);
y(aa) = d*(1.0-abs(x(aa)));

