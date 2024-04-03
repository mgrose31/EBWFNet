function [zz]=AnnularZernike(m, n, r, theta, e)
% SYNTAX:
% [zz]=AnnularZernike(m, n, r, theta, e)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% AnnularZernike(m, n, r, theta) = the (m,n)th circular Zernike polynomial 
% for radial coordinate r and azimuthal angle theta. The polynomial is 
% defined to be zero for r > 1. The circular Zernike polynomials form a 
% complete orthonormal set over the unit disc.
%
% AnnularZernike(m, n, r, theta, e) = the (m,n)th annular Zernike 
% polynomial, following the definition of Mahajan (JOSA 71, 75, 1981), for 
% obscuration e. The polynomial is defined to be zero for r > 1 or r < e. 
% The annular Zernike polynomials form a complete orthonormal set over the 
% annular domain e < r < 1.
% 
% Keith@MZA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REFERENCES:
% 1) V.N. Mahajan, J. Opt. Soc. Am., Vol 71 No 1, 75-85, Jan 1981; ...
%    Errata, J. Opt. Soc. Am., Vol A1, 685, 1984
% 2) V.N. Mahajan, J. Opt. Soc. Am. A., Vol 1 No 6, 685, June 1984
% 3) V.N. Mahajan, Suppliment to Applied Optics, Vol 5 No 11, 8125-8127, ...
%    Dec 1994
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: AnnularZernike.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

muse=abs(m);
zz=sqrt((n+1)/pi).*R(muse, n, r, e);
if(m<=0)
    zz=zz.*cos(m*theta);
else
    zz=zz.*sin(m*theta);
end
return
%
function [zz] = R(m, n, r, e)
% disp(['R(',num2str(m),', ',num2str(n),', ',num2str([r(1) r(end)]),', ',num2str(e),')'])
m=abs(m);
%
if(e==0)
    % Circular Zernike
    zz=0.0.*r;
    for s=0:(n-m)/2
        zz=zz+aGZ(m,n,s).*r.^(n-2*s);
    end
elseif(m==0 & mod(n,2)==0)
    % Ref 1, Equation (18); Ref 3, Equation (6)
    u=(r.^2 - e^2)/(1 - e^2);
    ind=find(u<0.0);u(ind)=0.0;
    zz = R(0, n, sqrt(u),0);
elseif(n==m)
    % Ref 1, Equation (25); Ref 3, Equation (7a)
    %     sig=0.0;
    %     for i=0:n
    %         sig=sig+e^(2*i);
    %     end
    %     zz = r.^n/sqrt(sig);
    % Ref 3, Equation (7b)
    zz = r.^n * sqrt((1 - e^2)/(1 - e^(2*n+2)));
elseif(m==n-2)
    % Ref 2, Equation (5); Ref 3, Equation (8)
    zz=( n*r.^n - (n-1)*( (1-e)^(2*n) / (1-e^(2*n-2)) ) * r.^(n-2) ) /...
        sqrt( (1-e^2)^(-1) * ( n^2 * (1-e^(2*n+2)) - (n^2-1)*(1-e^(2*n))^2 / (1-e^(2*n-2)) ) );
elseif(m==1 & n==5)
    % Ref 2, Equation (1); Ref 3, Table 1
    zz =( 10*(1 + 4*e^2 + e^4)*r.^5 - 12*(1 + 4*e^2 + 4*e^4 + e^6)*r.^3 + 3*(1 + 4*e^2 + 10*e^4 + 4*e^6 + e^8)*r ) /...
        ( (1 - e^2)^2 * sqrt((1 + 4*e^2 + e^4)*(1 + 9*e^2 + 9*e^4 + e^6)) );
elseif(m==2 & n==6)
    % Ref 2, Equation (4); Ref 3, Table 1
    zz =( 15*(1 + 4*e^2 + 10*e^4 + 4*e^6 + e^8)*r.^6 - 20*(1 + 4*e^2 + 10*e^4 + 10*e^6 + 4*e^8 + e^10)*r.^4 ...
        + 6*(1 + 4*e^2 + 10*e^4 + 20*e^6 + 10*e^8 + 4*e^10 + e^12)*r.^2 ) /...
        ( (1 - e^2 )^2 * sqrt((1 + 4*e^2 + 10*e^4 + 4*e^6 + e^8)*(1 + 9*e^2 + 45*e^4 + 65*e^6 + 45*e^8 + 9*e^10 + e^12)) );
elseif(m==3 & n==7)
    % ToMatlab[Simplify[RadialNM[7,3,r,e]]]
    zz =((-1)+e.^2).^(-2).*((1+4.*e.^2+10.*e.^4+20.*e.^6+10.*e.^8+4.* ...
        e.^10+e.^12).*(1+9.*e.^2+45.*e.^4+165.*e.^6+270.*e.^8+270.*e.^10+ ...
        165.*e.^12+45.*e.^14+9.*e.^16+e.^18)).^(-1/2).*r.^3.*(10+10.* ...
        e.^16+(-30).*r.^2+21.*r.^4+e.^14.*(40+(-30).*r.^2)+e.^12.*(100+( ...
        -120).*r.^2+21.*r.^4)+4.*e.^10.*(50+(-75).*r.^2+21.*r.^4)+10.* ...
        e.^8.*(35+(-60).*r.^2+21.*r.^4)+4.*e.^2.*(10+(-30).*r.^2+21.*r.^4) ...
        +10.*e.^4.*(10+(-30).*r.^2+21.*r.^4)+20.*e.^6.*(10+(-30).*r.^2+ ...
        21.*r.^4));
elseif(m==1 & n==7)
    % ToMatlab[Simplify[RadialNM[7,1,r,e]]]
    zz =((-1)+e.^2).^(-3).*((1+9.*e.^2+9.*e.^4+e.^6).*(1+16.*e.^2+36.* ...
        e.^4+16.*e.^6+e.^8)).^(-1/2).*r.*(4+4.*e.^12+(-30).*r.^2+60.*r.^4+ ...
        (-35).*r.^6+(-6).*e.^10.*((-6)+5.*r.^2)+30.*e.^8.*(6+(-9).*r.^2+ ...
        2.*r.^4)+(-5).*e.^6.*((-52)+150.*r.^2+(-108).*r.^4+7.*r.^6)+(-15) ...
        .*e.^4.*((-12)+50.*r.^2+(-60).*r.^4+21.*r.^6)+(-9).*e.^2.*((-4)+ ...
        30.*r.^2+(-60).*r.^4+35.*r.^6));
elseif(m==4 & n==8)
    % ToMatlab[Simplify[RadialNM[8,4,r,e]]]
    zz =((-1)+e.^2).^(-2).*((1+4.*e.^2+10.*e.^4+20.*e.^6+35.*e.^8+20.* ...
        e.^10+10.*e.^12+4.*e.^14+e.^16).*(1+9.*e.^2+45.*e.^4+165.*e.^6+ ...
        495.*e.^8+846.*e.^10+994.*e.^12+846.*e.^14+495.*e.^16+165.*e.^18+ ...
        45.*e.^20+9.*e.^22+e.^24)).^(-1/2).*r.^4.*(15+15.*e.^20+(-42).* ...
        r.^2+28.*r.^4+e.^18.*(60+(-42).*r.^2)+35.*e.^12.*(15+(-24).*r.^2+ ...
        8.*r.^4)+70.*e.^10.*(12+(-21).*r.^2+8.*r.^4)+2.*e.^16.*(75+(-84).* ...
        r.^2+14.*r.^4)+4.*e.^14.*(75+(-105).*r.^2+28.*r.^4)+4.*e.^2.*(15+( ...
        -42).*r.^2+28.*r.^4)+10.*e.^4.*(15+(-42).*r.^2+28.*r.^4)+20.* ...
        e.^6.*(15+(-42).*r.^2+28.*r.^4)+35.*e.^8.*(15+(-42).*r.^2+28.* ...
        r.^4));
elseif(m==2 & n==8)
    % ToMatlab[Simplify[RadialNM[8,2,r,e]]]
    zz =((-1)+e.^2).^(-3).*((1+9.*e.^2+45.*e.^4+65.*e.^6+45.*e.^8+9.* ...
        e.^10+e.^12).*(1+16.*e.^2+136.*e.^4+416.*e.^6+626.*e.^8+416.* ...
        e.^10+136.*e.^12+16.*e.^14+e.^16)).^(-1/2).*r.^2.*(10+10.*e.^18+( ...
        -60).*r.^2+105.*r.^4+(-56).*r.^6+e.^16.*(90+(-60).*r.^2)+15.* ...
        e.^14.*(30+(-36).*r.^2+7.*r.^4)+e.^10.*(2700+(-6900).*r.^2+4725.* ...
        r.^4+(-504).*r.^6)+e.^12.*(1650+(-2700).*r.^2+945.*r.^4+(-56).* ...
        r.^6)+(-9).*e.^2.*((-10)+60.*r.^2+(-105).*r.^4+56.*r.^6)+(-45).* ...
        e.^4.*((-10)+60.*r.^2+(-105).*r.^4+56.*r.^6)+(-15).*e.^8.*((-180)+ ...
        600.*r.^2+(-595).*r.^4+168.*r.^6)+(-5).*e.^6.*((-330)+1380.*r.^2+( ...
        -1785).*r.^4+728.*r.^6));      
else
    disp(['Error: Don"t know how to calculate radial polnomial R_',int2str(n),'^',int2str(m),'; returning zeros'])
    zz=0.0.*r;
end
ind=find(r<e | r>1);
zz(ind)=0.0;
% disp(['R(',num2str(m),', ',num2str(n),', ',num2str([r(1) r(end)]),', ',num2str(e),')= ',num2str([zz(1) zz(end)])])
return