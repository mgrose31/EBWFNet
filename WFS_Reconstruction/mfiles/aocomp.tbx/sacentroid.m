function [centroid, peakpix, maxval] = sacentroid(sa, relthresh, thresh, type)
% SYNTAX:
% [centroid, peakpix, maxval] = sacentroid(sa, relthresh, thresh, type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% sa [ ] = 
% relthresh [ ] = 
% thresh [ ] = 
% type [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% centroid [ ] = 
% peakpix [ ] = 
% maxval [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: sacentroid.m 3577 2016-08-26 17:04:29Z rpraus $

%% BEGIN_CODE

if (~exist('relthresh','var')), relthresh = 0; end;
if (isempty(relthresh)), relthresh = 0; end;
if (~exist('thresh','var')), thresh = 0; end;
if (isempty(thresh)), thresh = 0; end;
if (~exist('type','var')), type = 1; end;
if (isempty(type)), type = 1; end;
if (relthresh > 0)
   [peakpix, maxval] = findpeakpix(sa);
   thresh = max([thresh, relthresh * maxval]);
end;
if (thresh > 0)
   if (type == 2)
%       ind = find(sa < thresh);
%       if (~isempty(ind))
%          sa(ind) = 0;
%       end;
      sa(sa < thresh) = 0;
   elseif (type == 1)
      sa = sa - thresh;
%       ind = find(sa < 0);
%       if (~isempty(ind))
%          sa(ind) = 0;
%       end;
      sa(sa < 0) = 0;
   elseif (type == 3)
      sa = sa - thresh;
%       ind = find(sa < 0);
%       if (~isempty(ind))
%          sa(ind) = 0;
%       end;
      ind = sa <= 0;
      sa(ind) = 0;
%       ind = find(sa > 0);
%       if (~isempty(ind))
%          sa(ind) = 1;
%       end;
      sa(~ind) = 1;
   else
      error('Bad centroid type.')   
   end;
end;
centroid = saqcentroid(sa);
if ((nargout > 1) && (~exist('peakpix','var')))
   [peakpix, maxval] = findpeakpix(sa);
end;