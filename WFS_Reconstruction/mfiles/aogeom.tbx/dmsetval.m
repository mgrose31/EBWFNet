function dmsetval
% SYNTAX:
% dmsetval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: dmsetval.m 3028 2010-09-21 21:04:58Z amoran $

%% BEGIN_CODE

dmglobal;

%set all these before plotting

dm_act_num = str2num(get(guid_act_num,'string'));
for ic = 1:2
   dm_act_spacing(ic) = str2num(get(guid_act_spacing(ic),'string'));
   dm_act_offset(ic) = str2num(get(guid_act_offset(ic),'string'));
   dm_act_center_offset(ic) = str2num(get(guid_act_center_offset(ic),'string'));
end
if (dm_act_spacing(2) == 0)
   dm_act_spacing(2) = dm_act_spacing(1);
end   
dochange = 0;
for ic = 1:3
   for jc = 1:2
      dm_act_radii(ic,jc) = str2num(get(guid_act_radii(ic,jc),'string'));
   end
end
for ic = 1:3
   if (dm_act_radii(ic,1) >= dm_act_radii(ic,2))
       dm_act_radii(ic,1) = 0.0;
       dochange = 1;
   end;
end
for ic = 2:-1:1
   if (dm_act_radii(ic,2) > dm_act_radii(ic+1,2))
      dm_act_radii(ic,2) = dm_act_radii(ic+1,2);
      dochange = 1;
   end;
end
for ic = 1:2
   if (dm_act_radii(ic,1) < dm_act_radii(ic+1,1))
      dm_act_radii(ic,1) = dm_act_radii(ic+1,1);
      dochange = 1;
   end;
end
if (dochange == 1)
   for ic = 1:3
      for jc = 1:2
         set(guid_act_radii(ic,jc), 'string', num2str(dm_act_radii(ic,jc)));
      end
   end
end;

dm_slave_corners_three_ways = get(guid_slave_corners_three_ways, 'value');

wfs_act_num = str2num(get(guiw_act_num,'string'));
for ic = 1:2
   wfs_subaper_offset(ic) = ...
      str2num(get(guiw_subaper_offset(ic),'string'));
end   
wfs_subaper_space = str2num(get(guiw_subaper_space,'string'));
wfs_subaper_clip = str2num(get(guiw_subaper_clip,'string'));
wfs_subaper_anulus = str2num(get(guiw_subaper_anulus,'string'));
if (wfs_subaper_anulus > wfs_subaper_clip)
   wfs_subaper_anulus = 0.0;
   set(guiw_subaper_anulus, 'string', num2str(wfs_subaper_anulus));
end;
wfs_full_only = get(guiw_full_only,'value');
wfs_num_steps = str2num(get(guiw_num_steps,'string'));


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
op_disp_counts = strcmp(get(guio_disp_slaves,'Label'),...
				'Hide Counts');

for ic = 1:5
    if (get(guii_function(ic),'value') == 1)
       dm_influence = ic;
    end
end
if (dm_influence == 1)
   dm_clamp_rad = str2num(get(guii_infl_width(1),'string')); 
   dm_ninfl = str2num(get(guii_infl_width(2),'string'));
else
   dm_infl_width(1) = str2num(get(guii_infl_width(1),'string')); 
   dm_infl_width(2) = str2num(get(guii_infl_width(2),'string'));
end
