function [val]=a(l,n,s)
% SYNTAX:
% [val]=a(l,n,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: a.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE 

sign=(-1)^s;
m=(n-l)/2.0;
val=sign*factorial2(n-s);    
val=val/(factorial2(s)*factorial2(m-s)*factorial2(n-m-s));
return;