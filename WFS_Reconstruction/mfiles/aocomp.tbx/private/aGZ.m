function [val]=aGZ(l,m,s)
% SYNTAX:
% [val]=aGZ(l,m,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: aGZ.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

sign=(-1)^s;
val=sign*factorial(m-s)/(factorial(s)*factorial((m+l)/2-s)*factorial((m-l)/2-s));
return;