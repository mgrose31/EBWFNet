function irrVarSW = aowfsGetIrrVarSW(frameInfo,wl);
% SYNTAX: 
% irrVarSW = aowfsGetIrrVarSW(frameInfo,wl);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:  
% aowfsGetIrrVarSW does sliding window averaging for irradiance data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% frameInfo [ ] =  The data structure containing per frame data for the 
%                  data set
% wl [ ] =  The window length in number of frames for the sliding window 
%           average
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPTUS:
% irrVarSW [ ] =  The sliding window averages of subaperture irradiances.
%                 The number of such values is the (number of frames - 
%                 windowlength) + 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: aowfsGetIrrVarSW.m 3051 2010-10-01 20:33:26Z amoran $

%% BEGIN_CODE
 
nfr = length(frameInfo);
for ii=1:1:nfr-wl+1
  r=0;
  for jj=1:1:wl
    rr = [frameInfo(ii+jj-1).irr];
    r = rr + r;
  end
  r = r/wl;
  r = reshape(r,1,prod(size(r)));
  rm = mean(r);
  rn = r/rm;
  irrVarSW(ii) = (std(rn))^2;
end

