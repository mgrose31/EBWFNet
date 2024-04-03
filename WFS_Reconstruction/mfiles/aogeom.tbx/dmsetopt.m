function dmsetopt(num_opt);
% SYNTAX:
% dmsetopt(num_opt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% num_opt [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: dmsetopt.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE

dmglobal;
if(num_opt == 1)
op_disp_subaper = ~strcmp(get(guio_disp_subaper,'Label'),...
				'Hide Subaperture Radius');
if(op_disp_subaper == 0)
   set(guio_disp_subaper,'Label', 'Show Subaperture Radius');
else
   set(guio_disp_subaper,'Label', 'Hide Subaperture Radius');
end
end


if(num_opt == 2)
op_disp_master = ~strcmp(get(guio_disp_master,'Label'),...
				'Hide Master Actuator Radius');
if(op_disp_master == 0)
   set(guio_disp_master,'Label', 'Show Master Actuator Radius');
else
   set(guio_disp_master,'Label', 'Hide Master Actuator Radius');
end
end


if(num_opt == 3)
op_disp_slave = ~strcmp(get(guio_disp_slave,'Label'),...
				'Hide Slave Actuator Radius');
if(op_disp_slave == 0)
   set(guio_disp_slave,'Label', 'Show Slave Actuator Radius');
else
   set(guio_disp_slave,'Label', 'Hide Slave Actuator Radius');
end
end

if(num_opt == 4)
op_disp_inert = ~strcmp(get(guio_disp_inert,'Label'),...
				'Hide Inert Actuator Radius');
if(op_disp_inert == 0)
   set(guio_disp_inert,'Label', 'Show Inert Actuator Radius');
else
   set(guio_disp_inert,'Label', 'Hide Inert Actuator Radius');
end
end


if(num_opt == 5)
op_disp_clamp = ~strcmp(get(guio_disp_clamp,'Label'),...
				'Hide Clamping Radius');
if(op_disp_clamp == 0)
   set(guio_disp_clamp,'Label', 'Show Clamping Radius');
else
   set(guio_disp_clamp,'Label', 'Hide Clamping Radius');
end
end


if(num_opt == 6)
op_disp_axes = ~strcmp(get(guio_disp_axes,'Label'),...
				'Hide Axes');
if(op_disp_axes == 0)
   set(guio_disp_axes,'Label', 'Show Axes');
else
   set(guio_disp_axes,'Label', 'Hide Axes');
end
end


if(num_opt == 7)
op_disp_nums = ~strcmp(get(guio_disp_nums,'Label'),...
				'Hide Actuator Numbers');
if(op_disp_nums == 0)
   set(guio_disp_nums,'Label', 'Show Actuator Numbers');
else
   set(guio_disp_nums,'Label', 'Hide Actuator Numbers');
end
end


if(num_opt == 8)
op_disp_slaves = ~strcmp(get(guio_disp_slaves,'Label'),...
				'Hide Slave Relationships');
if(op_disp_slaves == 0)
   set(guio_disp_slaves,'Label', 'Show Slave Relationships');
else
   set(guio_disp_slaves,'Label', 'Hide Slave Relationships');
end
end


if(num_opt == 9)
op_disp_counts = ~strcmp(get(guio_disp_counts,'Label'),'Hide Counts');
if(op_disp_counts == 0)
   set(guio_disp_counts,'Label', 'Show Counts');
else
   set(guio_disp_counts,'Label', 'Hide Counts');
end
end


dmplot(0);
