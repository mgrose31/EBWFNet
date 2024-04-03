function [wg, dc] = satinker2d(frame, wgin, relthresh, thresh, type, mask)
% SYNTAX: 
% [wg, dc] = satinker2d(frame, wgin, relthresh, thresh, type, mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: satinker2d.m 3027 2010-09-21 21:04:10Z amoran $

%% BEGIN_CODE

if (~exist('relthresh','var')), relthresh = 0; end;
if (isempty(relthresh)), relthresh = 0; end;
if (~exist('thresh','var')), thresh = 0; end;
if (isempty(thresh)), thresh = 0; end;
if (~exist('type','var')), type = 1; end;
if (isempty(type)), type = 1; end;
if (~exist('mask','var')), mask = []; end;
if (isempty(mask)), mask = []; end;
dc = fmins('sarmsdisp', [0,0], [], [], frame, wgin, relthresh, thresh, type, [], mask);
wg = wpgtrans(wgin, dc(1), dc(2));
