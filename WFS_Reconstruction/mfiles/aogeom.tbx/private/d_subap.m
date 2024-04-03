function d_subap
% SYNTAX:
% d_subap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: d_subap.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE

dmglobal;

subaprad = 0.5*wfs_full_only*sqrt(wfs_subaper_space^2+wfs_subaper_space^2);
comp_nsubap = 0;
for xc = 1:wfs_act_num
   for yc = 1:wfs_act_num
      d = ((comp_mx_wfsx(yc,xc))^2 + (comp_mx_wfsy(yc,xc))^2);
      if ((d <= (wfs_subaper_clip - subaprad)^2) & ...
          (d >= (wfs_subaper_anulus + subaprad)^2))
         comp_nsubap = comp_nsubap+ 1;
         % This is to define a contour:
         xpoints = [1 -1 -1  1 1  1 -1 -1 1];
         ypoints = [1  1 -1 -1 1 -1 -1  1 1];
         xpoints = xpoints*wfs_subaper_space*0.4+ comp_mx_wfsx(yc,xc);
         ypoints = ypoints*wfs_subaper_space*0.4+ comp_mx_wfsy(yc,xc);
         tmp_h = patch('XData', xpoints,'YData', ypoints,...
                       'EdgeColor', 'red', 'FaceColor', 'red');
         % This is dark Red Color that I used for coloring subap
         % when it was represented as afilled square:  [0.44313, 0, 0.06666]);
      end
   end
end
