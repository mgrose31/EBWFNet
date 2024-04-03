function loadcomp
% SYNTAX:
% loadcomp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% The purpose of this function is to compute what is not yet defined
% right after loading a geometry file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: loadcomp.m 3028 2010-09-21 21:04:58Z amoran $

%% BEGIN_CODE

dmglobal;

comp_mx_acttype = zeros(dm_act_num,dm_act_num);
% master = 1, slave = 2, inert = 3, nonexistant = 0;
comp_mx_actx = zeros(dm_act_num,dm_act_num);
comp_mx_acty = zeros(dm_act_num,dm_act_num);
comp_mx_actnum = zeros(dm_act_num,dm_act_num);
comp_vec_actnum = zeros(dm_act_num^2,2);

comp_nact = 0;
comp_nslav = 0;
comp_nmas = 0;

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

        if (comp_mx_acttype(yc,xc) == 1)
           comp_nmas = comp_nmas + 1;
        elseif (comp_mx_acttype(yc,xc) == 2)
           comp_nslav = comp_nslav + 1;
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
