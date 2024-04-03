function dm_redisplay_options;
% SYNTAX:
% dm_redisplay_options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: dm_redisplay_options.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE

op_disp_subaper = strcmp(get(guio_disp_subaper,'Label'),...
				'Hide Subaperature Radius')
op_disp_master = strcmp(get(guio_disp_master,'Label'),...
				'Hide Master Actuator Radius');
op_disp_slave = strcmp(get(guio_disp_slave,'Label'),...
				'Hide Slave Actuator Radius');
op_disp_inert = strcmp(get(guio_disp_inert,'Label'),...
				'Hide Inert Actuator Radius');
op_disp_clamp = strcmp(get(guio_disp_clamp,'Label'),...
				'Hide Clamping Radius');
op_disp_axes = strcmp(get(guio_disp_axes,'Label'),...
				'Hide Axes');
op_disp_nums = strcmp(get(guio_disp_nums,'Label'),...
				'Hide Actuator Numbers');
op_disp_slaves = strcmp(get(guio_disp_slaves,'Label'),...
				'Hide Slave Relationships');
dmplot(0);
