function dmreset
% SYNTAX:
% dmreset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: dmreset.m 3028 2010-09-21 21:04:58Z amoran $

%% BEGIN_CODE

dmglobal;

%reset gui values to true values

set(guid_act_num,'string',num2str(dm_act_num)); 
for ic = 1:2
   set(guid_act_spacing(ic),'string',num2str(dm_act_spacing(ic))); 
   set(guid_act_offset(ic),'string',num2str(dm_act_offset(ic))); 
   set(guid_act_center_offset(ic),'string',num2str(dm_act_center_offset(ic))); 
end
   
for ic = 1:3
   for jc = 1:2
      set(guid_act_radii(ic,jc),'string',num2str(dm_act_radii(ic,jc))); 
   end
end
set(guid_slave_corners_three_ways, 'value', dm_slave_corners_three_ways);


set(guiw_act_num,'string',num2str(wfs_act_num)); 
for ic = 1:2
   set( guiw_subaper_offset(ic),'string',...
	num2str(wfs_subaper_offset(ic)));
end   
set(guiw_subaper_space,'string',num2str(wfs_subaper_space)); 
set(guiw_subaper_clip,'string',num2str(wfs_subaper_clip)); 
set(guiw_subaper_anulus,'string',num2str(wfs_subaper_anulus)); 
set(guiw_full_only,'value',wfs_full_only); 
set(guiw_num_steps,'string',num2str(wfs_num_steps)); 


if(op_disp_subaper == 1)
  set(guio_disp_subaper, 'Label', 'Hide Subaperture Radius');
else
  set(guio_disp_subaper, 'Label', 'Show Subaperture Radius');
end

if(op_disp_master == 1)
  set(guio_disp_master, 'Label', 'Hide Master Actuator Radius');
else
  set(guio_disp_master, 'Label', 'Show Master Actuator Radius');
end

if(op_disp_slave == 1)
  set(guio_disp_slave, 'Label', 'Hide Slave Actuator Radius');
else
  set(guio_disp_slave, 'Label', 'Show Slave Actuator Radius');
end

if(op_disp_inert == 1)
  set(guio_disp_inert, 'Label', 'Hide Inert Actuator Radius');
else
  set(guio_disp_inert, 'Label', 'Show Inert Actuator Radius');
end

if(op_disp_clamp == 1)
  set(guio_disp_clamp, 'Label', 'Hide Clamping Radius');
else
  set(guio_disp_clamp, 'Label', 'Show Clamping Radius');
end

if(op_disp_axes == 1)
  set(guio_disp_axes, 'Label', 'Hide Axes');
else
  set(guio_disp_axes, 'Label', 'Show Axes');
end

if(op_disp_nums == 1)
  set(guio_disp_nums, 'Label', 'Hide Actuator Numbers');
else
  set(guio_disp_nums, 'Label', 'Show Actuator Numbers');
end

if(op_disp_slaves == 1)
  set(guio_disp_slaves, 'Label', 'Hide Slave Relationships');
else
  set(guio_disp_slaves, 'Label', 'Show Slave Relationships');
end

if(op_disp_counts == 1)
  set(guio_disp_counts, 'Label', 'Hide Counts');
else
  set(guio_disp_counts, 'Label', 'Show Counts');
end

for ic = 1:5
   if (dm_influence == ic)
      set(guii_function(ic),'value',1);
   else
      set(guii_function(ic),'value',0);
   end
end
if (dm_influence == 1)
   set(guii_infl_width(1),'string',num2str(dm_clamp_rad)); 
   set(guii_infl_width(2),'string',num2str(dm_ninfl));
else
   set(guii_infl_width(1),'string',num2str(dm_infl_width(1))); 
   set(guii_infl_width(2),'string',num2str(dm_infl_width(2)));
end

tempstr = get(guii_function(dm_influence),'userdata');

set(guii_infl_label(1),'string', tempstr(1,:));
set(guii_infl_label(2),'string', tempstr(2,:));

