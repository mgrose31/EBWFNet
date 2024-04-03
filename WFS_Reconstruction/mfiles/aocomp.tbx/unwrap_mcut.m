function [surf] = unwrap_mcut(phase,mask,qual,mode,tsize,surf_fig)
% SYNTAX:
% [surf] = unwrap_mcut(phase,mask,qual,mode,tsize,surf_fig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Mask cut algorithm for unwrapping two-dimensional phase.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WHERE:
%   phase is the giving wrapped phase matrix, must be real.
%   mask is the mask matrix, or [].  qual is the quality value
%   or correlation data, it can be either a matrix or [].
%   mode is mkey' is a keyword designating the source of the quality values
%   that guide the masking path.  The value of mkey must be:
%   'none'      or    0  for the default mode;
%   'min_grad'  or    1  for minimum gradients;
%   'min_var'   or    2  for minimum phase variances;
%   'max_pseu'  or    3 for maximum pseudocorrelations (default is none); 
%   'max_corr'  or    4 for maximum cross correlations (i.e., corr data).
%   tsize is the size of the square template for averaging the corr data
%   or quality values, the default value is 1. It must be either a
%   number or [].
% Tempus: 
% Grid<float> = unwrap_mcut(Grid<float>, Grid<float>, Grid<float>, int, int, int)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% phase [ ] =
% mask [ ] = 
% qual [ ] = 
% mode [ ] = 
% tsize [ ] =
% surf_fig [ ] =
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPTUS:
% surf [ ] =
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: unwrap_mcut.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

mkey = [];
if (isempty(mode) ~= 1)
	if (ischar(mode) == 1)
   	switch lower(mode)
   	case 'none'
      	mkey = 0;
   	case 'min_grad'
      	mkey = 1;
   	case 'min_var'
      	mkey = 2;
   	case 'max_pseu'
      	mkey = 3;
   	case 'max_corr'
      	mkey = 4;
   	otherwise
      	disp('Unknown method.');
      	mkey = 0;
   	end
	else
   	mkey = mode;
   end
end

if (~isempty(mask))
   the_mask = mask.g;
else
  ii = find(phase.g);
  the_mask = zeros(size(phase.g));
  the_mask(ii) = 1;
end

if (~isempty(qual))
   the_qual = qual.g;
end

if (isempty(tsize))
   tsize = 1;
end

surface = mexmcut(phase.g,the_mask,the_qual,mkey,tsize);
surf.x = [1:size(surface,1)];
surf.y = [1:size(surface,2)];
surf.g = surface;
figure(surf_fig);
imagesc(surf.x, surf.y, surf.g);
axis image;
colorbar;
drawnow;
