function [surf]=unwrap_gold(phase,mask,surf_fig)
% SYNTAX:
% [surf]=unwrap_gold(phase,mask,surf_fig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% GOLD Goldstein's algorithm for unwrapping two-dimensional phase.
% Where:
%   phase is the giving wrapped phase matrix, must be real.
%   mask is the mask matrix, it's real too. Also mask and phase
%   must be the same in size.  If there is no mask, use [].
% Tempus: Grid<float> = unwrap_gold(Grid<float>, Grid<float>, int)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% phase [ ] =
% mask [ ] = 
% surf_fig [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% surf [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: unwrap_gold.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

if (~isempty(mask))
  the_mask = mask.g;
else
  ii = find(phase.g);
  the_mask = zeros(size(phase.g));
  the_mask(ii) = 1;
end

the_surf = mexgold(phase.g,the_mask);
surf.x = [1:size(the_surf,1)];
surf.y = [1:size(the_surf,2)];
surf.g = the_surf;
figure(surf_fig);
imagesc(surf.x, surf.y, surf.g);
axis image;
colorbar;
drawnow;
