function dmdefaul
% SYNTAX:
% dmdefaul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: dmdefaul.m 3028 2010-09-21 21:04:58Z amoran $

%% BEGIN_CODE

dmglobal;

dm_act_num = 7;
dm_act_spacing = [0.1,0.1];
dm_act_offset = [0,0];
dm_act_center_offset = [0,0];

% this is master(inner_r, outer_r), slave(inner_r, outer_r), 
% inert(inner_r, outer_r,)
dm_act_radii = [0.2 0.35; 0.0 0.4 ; 0.0 0.45];
dm_slave_corners_three_ways = 1;

wfs_act_num = 8;
wfs_subaper_offset = [0.0,0.0];
wfs_subaper_space = 0.1;
wfs_subaper_clip = 0.28;
wfs_subaper_anulus = 0;
wfs_full_only = 0;
wfs_num_steps = 6;

op_filename = 'dmgeom';
op_disp_subaper = 1;
op_disp_master = 1;
op_disp_slave = 1;
op_disp_inert = 1;
op_disp_clamp = 1;
op_disp_axes = 1;
op_disp_nums = 0;
op_disp_slaves = 1;
op_disp_counts = 0;

dm_influence = 1;
dm_clamp_rad = 0.5;
dm_ninfl = 2;
dm_infl_width = dm_act_spacing;

%Colors used in GUI:
%BGColor = [.5, .5, .5]; % black
%labelFGColor = [.0, .0, .0]; % some kind of gray
BGColor = [.0, .0, .0]; % black
labelFGColor = [.7, .7, .7]; % some kind of gray

SlaveColor = 'magenta';
MasterColor = 'white';
InertColor  = 'green';


 h_inertCircles = -1;
 h_clampCircle = -1;
 h_slaveCircles = -1;
 h_masterCircles = -1;
 h_subapClipCircle = -1;
 h_subapAnulusCircle = -1;
 h_slave2masterLines = -1;
 h_actuators = -1;

%Things that we need to keep for the GUI only purposes:

 lastSelectedObject = [ -1 -1 -1 -1];
 isMouseButtonDown = 0;
 lastClickTime = -1;
 doubleClickInterval = 0.3;


%Some constants:
 SUPPORTS_DOUBLE_CLICK = 1;

%Types of objects that can be selected:
INERT_INNER_CIRCLE = 1;
INERT_OUTER_CIRCLE = 2;

SLAVE_INNER_CIRCLE = 3;
SLAVE_OUTER_CIRCLE = 4;

MASTER_INNER_CIRCLE = 5;
MASTER_OUTER_CIRCLE = 6;

CLAMPING_CIRCLE = 7;

SUBAPERTURE_INNER_CIRCLE = 8;
SUBAPERTURE_OUTER_CIRCLE = 9;

SLAVE_TO_MASTER_LINE = 10;

ACTUATOR = 11;
