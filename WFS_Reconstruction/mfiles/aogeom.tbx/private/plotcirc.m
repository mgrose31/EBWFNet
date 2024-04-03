function plotcirc
% SYNTAX:
% plotcirc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: plotcirc.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE

dmglobal;

if (op_disp_clamp & (dm_influence == 1))
   h_clampCircle = d_circle(dm_clamp_rad,dm_act_offset(1),dm_act_offset(2),'cyan');
end
if (op_disp_inert)
   h_inertCircles = dmcircle(dm_act_radii(3,:),dm_act_center_offset(1),dm_act_center_offset(2),'green');
end
if (op_disp_slave)
   h_slaveCircles = dmcircle(dm_act_radii(2,:),dm_act_center_offset(1),dm_act_center_offset(2),'magenta');
end
if (op_disp_master)
   h_masterCircles = dmcircle(dm_act_radii(1,:),dm_act_center_offset(1),dm_act_center_offset(2),'white');
end
if (op_disp_subaper) 
   h_subapClipCircle = d_circle(wfs_subaper_clip, 0, 0,'red');
%   h_subapClipCircle = d_circle(wfs_subaper_clip, wfs_subaper_offset(1), wfs_subaper_offset(2),'red');
   if(wfs_subaper_anulus > 0)
      h_subapAnulusCircle = d_circle(wfs_subaper_anulus, 0, 0, 'red');
%      h_subapAnulusCircle = d_circle(wfs_subaper_anulus, wfs_subaper_offset(1), wfs_subaper_offset(2), 'red');
   end
end
