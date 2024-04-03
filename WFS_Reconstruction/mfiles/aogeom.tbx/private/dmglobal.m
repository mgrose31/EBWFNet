global dm_act_num;
global dm_act_spacing;
global dm_act_offset;
global dm_act_center_offset;

global dm_act_radii;
global dm_slave_corners_three_ways;

global wfs_act_num;
global wfs_subaper_offset;
global wfs_subaper_space;
global wfs_subaper_clip;
global wfs_subaper_anulus;
global wfs_full_only;
global wfs_num_steps;

global op_filename;
global op_disp_subaper;
global op_disp_master;
global op_disp_slave;
global op_disp_inert;
global op_disp_clamp;
global op_disp_axes;
global op_disp_nums;
global op_disp_slaves;
global op_disp_counts;

global gui_rectangle;
global gui_main;
global gui_box;

global guid_act_num;
global guid_act_spacing;
global guid_act_offset;
global guid_act_center_offset;
global guid_act_radii;
global guid_slave_corners_three_ways;

global guiw_act_num;
global guiw_subaper_offset;
global guiw_subaper_space;
global guiw_subaper_clip;
global guiw_subaper_anulus;
global guiw_full_only;
global guiw_num_steps;

global guio_filename;
global guio_disp_subaper;
global guio_disp_master;
global guio_disp_slave;
global guio_disp_inert;
global guio_disp_clamp;
global guio_disp_axes;
global guio_disp_nums;
global guio_disp_slaves;
global guio_disp_counts;

global dm_influence;
global dm_infl_width;
global dm_clamp_rad;
global dm_ninfl;

global guii_function;
global guii_infl_width;
global guii_infl_label;

global comp_nact;
global comp_nslav;
global comp_nmas;
global comp_nmsmax;
global comp_mx_acttype;
global comp_mx_actx;
global comp_mx_acty;
global comp_mx_actnum;
global comp_vec_actnum;
global comp_mx_wfsx;
global comp_mx_wfsy;
global comp_slave_defs;
global comp_slave_weights;
global comp_nsubap;

global BGColor;
global labelFGColor;

global SlaveColor;
global MasterColor;
global InertColor;


%Object handles:
global h_inertCircles;
global h_clampCircle;
global h_slaveCircles;
global h_masterCircles;
global h_subapClipCircle;
global h_subapAnulusCircle;

global h_slave2masterLines; % holds slave actuators with all 
			    % handles of lines that go from them
global h_Lines;		%Holds handles of lines and numbers of actuators it connects.

global h_actuators;

%Things that we need to keep for the GUI only purposes:

global lastSelectedObject;
global isMouseButtonDown;
global lastClickTime;
global doubleClickInterval;

%Some constants:
global SUPPORTS_DOUBLE_CLICK;

%Types of objects that can be selected:
global INERT_INNER_CIRCLE;
global INERT_OUTER_CIRCLE;

global SLAVE_INNER_CIRCLE;
global SLAVE_OUTER_CIRCLE;

global MASTER_INNER_CIRCLE;
global MASTER_OUTER_CIRCLE;

global CLAMPING_CIRCLE;

global SUBAPERTURE_INNER_CIRCLE;
global SUBAPERTURE_OUTER_CIRCLE;

global SLAVE_TO_MASTER_LINE;

global ACTUATOR;
