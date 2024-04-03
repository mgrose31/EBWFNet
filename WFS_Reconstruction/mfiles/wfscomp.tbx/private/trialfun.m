function f=trialfun(x,d1,d2);
% SYNTAX:
% f=trialfun(x,d1,d2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: trialfun.m 3051 2010-10-01 20:33:26Z amoran $

%% BEGIN_CODE

f=sum((d1-x*d2).^2);

