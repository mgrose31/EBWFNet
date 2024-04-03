function [imin, imax, simg, ncol, imghndl] = tvphs(c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: tvphs.m 3051 2010-10-01 20:33:26Z amoran $

%% BEGIN_CODE

[imin, imax, simg, ncol, imghndl] = tvscl(atan2(imag(c),real(c)));
