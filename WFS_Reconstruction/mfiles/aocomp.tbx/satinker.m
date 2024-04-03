function [wg, dc] = satinker(frame, wgin, relthresh, thresh, type, mask, method)
% SYNTAX: 
% [wg, dc] = satinker(frame, wgin, relthresh, thresh, type, mask, method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% frame [ ] = 
% wgin [ ] = 
% relthresh [ ] = 
% thresh [ ] = 
% type [ ] = 
% mask [ ] = 
% method [ ] =
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% wg [ ] = 
% dc [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: satinker.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

if (~exist('relthresh','var')), relthresh = 0; end;
if (isempty(relthresh)), relthresh = 0; end;
if (~exist('thresh','var')), thresh = 0; end;
if (isempty(thresh)), thresh = 0; end;
if (~exist('type','var')), type = 1; end;
if (isempty(type)), type = 1; end;
if (~exist('mask','var')), mask = []; end;
if (isempty(mask)), mask = []; end;
if (~exist('method','var')), method = 3; end;
if (isempty(method)), method = 3; end;
if ((method == 1) | (method == 3))
   [wg1, dc1] = satinker1d(frame, wgin, relthresh, thresh, type, mask);
else
   wg1 = wgin;
end;
if (method > 1)
   [wg2, dc2] = satinker2d(frame, wg1, relthresh, thresh, type, mask);
end;
if (exist('dc1','var'))
   if (exist('dc2','var'))
      dc = dc1 + dc2;
      wg = wpgtrans(wgin, dc(1), dc(2));
   else
      dc = dc1;
      wg = wg1;
   end;
else
   dc = dc2;
   wg = wg2;
end;
