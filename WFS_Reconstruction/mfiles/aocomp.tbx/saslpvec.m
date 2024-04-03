function sv = saslpvec(pos, mask)
% SYNTAX:
% sv = saslpvec(pos, mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS: 
% pos [ ] = 
% mask [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% sv [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: saslpvec.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

ind = find(mask > 0);
n = prod(size(ind));
sv = zeros(n,2);
slps = squeeze(pos(2,:,:));
sv(:,1) = slps(ind);
slps = squeeze(pos(1,:,:));
sv(:,2) = slps(ind);
sv = reshape(sv, 2*n, 1);

