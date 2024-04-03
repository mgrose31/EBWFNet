function [pp, mv] = findpeakpix(sa)
% SYNTAX:
% [pp, mv] = findpeakpix(sa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% sa [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS: 
% pp [ ] = 
% mv [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: findpeakpix.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

[sa1 sa2 sa3] = size(sa);
sa = reshape(sa,sa1*sa2,sa3);
[mv in] = max(sa,[],1);
[in1 in2] = ind2sub([sa1,sa2],in);
pp = [in1(:), in2(:)];
mv = mv(:);
