function [n,m,sincos]=pzrntonmv(k);
% SYNTAX:
% [n,m,sincos]=pzrntonmv(k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Computes a mapping from an indexing of the Zernike polynomials by Z to 
% the double indexing of the radial polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: pzrntonmv.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE
 
sincos=ones(1,length(k));
n=fix(sqrt(2*k+.25)-.5);
aa = (n.*(n-1))/2;
iq = mod(aa,2);
m=2*fix((k+1-iq)/2)-(n.*(n+1))/2+iq;
isn = mod(k-aa,2);
xx =find((m==0) | (isn == 1));
sincos(xx) = 0;

