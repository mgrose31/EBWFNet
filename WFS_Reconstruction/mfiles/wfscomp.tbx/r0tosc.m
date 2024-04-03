function [rout,sig,ssout]=r0tosc(s,procInfo);
% SYNTAX:
% [rout,sig,ssout]=r0tosc(s,procInfo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% analyze the r0 of the data contained in outz 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: r0tosc.m 3051 2010-10-01 20:33:26Z amoran $

%% BEGIN_CODE
 
nsimp = 6;
iter=5;
guess=10.;
nz = size(procInfo.w,2);
[w] = mkwv(procInfo,nz,nsimp);
radap = procInfo.outradius;
[c] = mknollv(nz,radap);
[rout,sig,ssout] = r0test(s,radap,w,c,guess,iter);
