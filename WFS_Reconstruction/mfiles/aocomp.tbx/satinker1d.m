function [wg, dc] = satinker1d(frame, wgin, relthresh, thresh, type, mask)
% SYNTAX: 
% [wg, dc] = satinker1d(frame, wgin, relthresh, thresh, type, mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: satinker1d.m 3027 2010-09-21 21:04:10Z amoran $

%% BEGIN_CODE

if (~exist('relthresh','var')), relthresh = 0; end;
if (isempty(relthresh)), relthresh = 0; end;
if (~exist('thresh','var')), thresh = 0; end;
if (isempty(thresh)), thresh = 0; end;
if (~exist('type','var')), type = 1; end;
if (isempty(type)), type = 1; end;
if (~exist('mask','var')), mask = []; end;
dcx = fmin('sarmsdisp', -max(wgin.pxb), max(wgin.pxf), [], frame, wgin, relthresh, thresh, type, 1, mask);
dc = [dcx, 0];
wg = wpgtrans(wgin, dc(1), dc(2));
dcy = fmin('sarmsdisp', -max(wg.pyb), max(wg.pyf), [], frame, wg, relthresh, thresh, type, 2, mask);
dc = [dcx, dcy];
wg = wpgtrans(wgin, dc(1), dc(2));
