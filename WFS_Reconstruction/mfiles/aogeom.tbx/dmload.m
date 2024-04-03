function dmload(filename)
% SYNTAX:
% dmload(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% filename [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% handle [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: dmload.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE

dmglobal;
%dmdefaul;

load(filename);

if (~exist('xact') | ~exist('yact') | ~exist('actype') |...
    ~exist('actnum') | ~exist('idxmas'))
   errordlg({'The specified file was not created by',...
             'AOGeom. The load will not be completed'},...
             'AOGeom: Error loading',...
             'modal');
   return;
end;   

if (exist('xact') & exist('yact') & exist('actype'))
   if (~exist('rmas_out'))
      rmas_out = 1.001 * sqrt(max((xact.^2 + yact.^2) .* (actype == 1)));
   end
   if (~exist('rslav_out'))
      rslav_out = 1.001 * sqrt(max((xact.^2 + yact.^2) .* (actype == -1)));
   end   
   if (~exist('ract_out'))
      ract_out = 1.001 * sqrt(max(xact.^2 + yact.^2));
   end

   if (~exist('rmas_in'))
      rmas_in = 0.999 * sqrt(min((xact.^2 + yact.^2) .* (actype == 1)));
      if(rmas_in < 0), rmas_in = 0; end
   end
   if (~exist('rslav_in'))
      rslav_in = 0.999 * sqrt(min((xact.^2 + yact.^2) .* (actype == -1)));
      if(rslav_in < 0), rslav_in = 0; end
   end   
   if (~exist('ract_in'))
      ract_in = 0.999 * sqrt(min(xact.^2 + yact.^2));
      if(ract_in < 0), ract_in = 0; end
   end
end

if (exist('rmas_out')),   dm_act_radii(1,2) = rmas_out; end
if (exist('rslav_out')),  dm_act_radii(2,2) = rslav_out; end
if (exist('ract_out')),   dm_act_radii(3,2) = ract_out; end

if (exist('rmas_in')),   dm_act_radii(1,1) = rmas_in; end
if (exist('rslav_in')),  dm_act_radii(2,1) = rslav_in; end
if (exist('ract_in')),   dm_act_radii(3,1) = ract_in; end


dm_act_offset(1:2) = 0.0;
if exist('xoff'), dm_act_offset(1) = xoff; end 
if exist('yoff'), dm_act_offset(2) = yoff; end 
dm_act_center_offset(1:2) = 0.0;
if exist('center_xoff'), dm_act_center_offset(1) = center_xoff; end 
if exist('center_yoff'), dm_act_center_offset(2) = center_yoff; end 
dm_slave_corners_three_ways = 0;
if exist('slave_corners_3ways')
   dm_slave_corners_three_ways = slave_corners_3ways;
end;

if exist('xact')
   for ic = 1:(nact-1)
       dm_act_spacing(1) = xact(ic+1) - xact(ic);
       if (dm_act_spacing(1) > 0)
          break;
       end
   end
end
if exist('yact')
   for ic = 1:(nact-1)
       dm_act_spacing(2) = yact(ic+1) - yact(ic);
       if (dm_act_spacing(2) > 0)
          break;
       end
   end
end
if exist('actnum')
   dm_act_num = actnum; 
else
   if (exist('xact'))
      dm_act_num = (max(xact) - min(xact))/dm_act_spacing(1) + 1;
      dm_act_num = round(dm_act_num);
   end
end

if exist('rsub_out')
   wfs_subaper_clip = rsub_out;
else
   if (exist('xsub') & exist('ysub'))
      rsub_out = 1.001 * sqrt(max(xsub.^2 + ysub.^2));
      wfs_subaper_clip = rsub_out;
   end
end
  
if exist('rsub_in')
   wfs_subaper_anulus = rsub_in;
else
   if (exist('xsub') & exist('ysub'))
      rsub_in = 0.999 * sqrt(min(xsub.^2 + ysub.^2));
      wfs_subaper_anulus = rsub_in;
   end
end
    
if exist('hs'), wfs_subaper_space = hs; end 
if exist('nsi'), wfs_num_steps = nsi; end 

if exist('xoffsub'), wfs_subaper_offset(1) = xoffsub; end 
if exist('yoffsub'), wfs_subaper_offset(2) = yoffsub; end 
if exist('fullsub'), wfs_full_only = fullsub; end
if exist('subnum')
   wfs_act_num = subnum; 
else
   if (exist('xsub'))
      wfs_act_num = (max(xsub) - min(xsub))/wfs_subaper_space + 1;
      wfs_act_num = round(wfs_act_num);
   end
end

if exist('dminfl');
   if (dminfl == 0); dm_influence = 1; end;
   if (dminfl == 4); dm_influence = 2; end;
   if (dminfl == 1); dm_influence = 3; end;
   if (dminfl == 2); dm_influence = 4; end;
   if (dminfl == 3); dm_influence = 5; end;
end;
if exist('rc'), dm_clamp_rad = rc; end 
if exist('ninfl'), dm_ninfl = ninfl; end
dm_infl_width = dm_act_spacing;
if exist('xwidth') dm_infl_width(1) = xwidth; end 
if exist('ywidth'), dm_infl_width(2) = ywidth; end 


comp_slave_defs = idxmas;
comp_slave_weights = weight;


comp_mx_acttype = zeros(dm_act_num,dm_act_num);
% master = 1, slave = 2, inert = 3, nonexistant = 0;
comp_mx_actx = zeros(dm_act_num,dm_act_num);
comp_mx_acty = zeros(dm_act_num,dm_act_num);
comp_mx_actnum = zeros(dm_act_num,dm_act_num);
comp_vec_actnum = zeros(dm_act_num^2,2);

comp_nact  = 0;
comp_nslav = nslav;
comp_nmas  = nmas;

first_coord = ((dm_act_num - 1) .* dm_act_spacing) ./ 2 .* (-1)...
              + dm_act_offset;
ic = 0;

for yc = 1:dm_act_num
    for xc = 1:dm_act_num
        comp_mx_actx(yc,xc) = first_coord(1) + dm_act_spacing(1)*(xc-1);
        comp_mx_acty(yc,xc) = first_coord(2) + dm_act_spacing(2)*(yc-1);

	distToCenter = ((comp_mx_actx(yc,xc) - dm_act_center_offset(1))^2 + ...
			(comp_mx_acty(yc,xc) - dm_act_center_offset(2))^2);

	% First we check if it is an nonexistant actuator
	if((distToCenter > (dm_act_radii(3,2)^2)) |...
	   (distToCenter < (dm_act_radii(3,1)^2)))
	  comp_mx_acttype(yc,xc) = 0;	
	else  %Otherwise we use the existing vector 'actype'
	  ic = ic + 1;
	  if(actype(ic) == -1)
	    comp_mx_acttype(yc,xc) = 2; % Slave
	  elseif(actype(ic) == 0)
	    comp_mx_acttype(yc,xc) = 3; % Inert
	  elseif(actype(ic) == 1)
	    comp_mx_acttype(yc,xc) = 1; % Master
	  end
	end	

        if (comp_mx_acttype(yc,xc) > 0)
           comp_nact = comp_nact + 1;
           comp_mx_actnum(yc,xc) = comp_nact;
           comp_vec_actnum(comp_nact,1) = xc;
           comp_vec_actnum(comp_nact,2) = yc;
        else
           comp_mx_actnum(yc,xc) = 0;
        end
    end
end

comp_mx_wfsx = zeros(wfs_act_num,wfs_act_num);
comp_mx_wfsy = zeros(wfs_act_num,wfs_act_num);
comp_nwfs = 0;

first_coord = ((wfs_act_num - 1) .* wfs_subaper_space) ./ 2 .* (-1);

for xc = 1:wfs_act_num
    for yc = 1:wfs_act_num
        comp_mx_wfsx(yc,xc) = first_coord + wfs_subaper_offset(1) + ...
               wfs_subaper_space*(xc-1);
        comp_mx_wfsy(yc,xc) = first_coord + wfs_subaper_offset(2) + ...
               wfs_subaper_space*(yc-1);
    end
end

dmreset;

s_tmp = size(comp_slave_defs);

comp_nmsmax = s_tmp(2);

