function rmsdisp = sarmsdisp(dc, frame, wgin, relthresh, thresh, type, ...
                   direction, mask)
% SYNTAX:
% rmsdisp = sarmsdisp(dc, frame, wgin, relthresh, thresh, type, ...
%           direction, mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS: 
% dc [ ] = 
% frame [ ] = 
% wgin [ ] = 
% relthresh [ ] = 
% thresh [ ] = 
% type [ ] = 
% direction [ ] = 
% mask [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% rmsdisp [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: sarmsdisp.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

if (~exist('relthresh','var')), relthresh = 0; end;
if (isempty(relthresh)), relthresh = 0; end;
if (~exist('thresh','var')), thresh = 0; end;
if (isempty(thresh)), thresh = 0; end;
if (~exist('type','var')), type = 1; end;
if (isempty(type)), type = 1; end;
if (prod(size(dc)) == 1)
   if (~exist('direction','var')), direction = 1; end;
   if (isempty(direction)), direction = 1; end;
   if (direction == 1)
      dc = [dc(1), 0];
   else
      dc = [0, dc(1)];
   end;
end;
if (~exist('mask','var')), mask = []; end;
[centroid, peakpix] = wfscentroid(frame, wpgtrans(wgin, dc(1), dc(2)), relthresh, thresh, type, mask);
if (type == 0)
   rmsdisp = sqrt(sum(sum(sum(peakpix.^2))));
else
   rmsdisp = sqrt(sum(sum(sum(centroid.^2))));
end;
