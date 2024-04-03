function nact = dmcomp;
% SYNTAX:
% nact = dmcomp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: dmcomp.m 3028 2010-09-21 21:04:58Z amoran $

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

first_coord = ((dm_act_num - 1) .* dm_act_spacing) ./ 2 .* (-1) + dm_act_offset;

for yc = 1:dm_act_num
   for xc = 1:dm_act_num
      comp_mx_actx(yc,xc) = first_coord(1) + dm_act_spacing(1)*(xc-1);
      comp_mx_acty(yc,xc) = first_coord(2) + dm_act_spacing(2)*(yc-1);

      % now we need to set (for all actuators) their type, depending on the 
      % location
      % First we compute distance^2 frmo the actuator to each of the centers
      % which could be the same point, but may differ.
   
      distToCenter = ((comp_mx_actx(yc,xc) - dm_act_center_offset(1))^2 + ...
                      (comp_mx_acty(yc,xc) - dm_act_center_offset(2))^2);

      % First we check if it is an nonexistant actuator
      if ((distToCenter > (dm_act_radii(3,2)^2)) | (distToCenter < (dm_act_radii(3,1)^2)))
         comp_mx_acttype(yc,xc) = 0;            
      elseif (((distToCenter >= (dm_act_radii(3,1)^2)) & (distToCenter < (dm_act_radii(2,1)^2))) |...
              ((distToCenter >= (dm_act_radii(2,2)^2)) & (distToCenter < (dm_act_radii(3,2)^2))))
         comp_mx_acttype(yc,xc) = 3; %inert
      elseif (((distToCenter >= (dm_act_radii(2,1)^2)) & (distToCenter < (dm_act_radii(1,1)^2))) |...
              ((distToCenter >= (dm_act_radii(1,2)^2)) & (distToCenter < (dm_act_radii(2,2)^2))))
         comp_mx_acttype(yc,xc) = 2; %slave
      elseif ((distToCenter >= (dm_act_radii(1,1)^2)) & (distToCenter < (dm_act_radii(1,2)^2)))
         comp_mx_acttype(yc,xc) = 1; %master
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
      comp_mx_wfsx(yc,xc) = first_coord + wfs_subaper_offset(1) + wfs_subaper_space*(xc-1);
      comp_mx_wfsy(yc,xc) = first_coord + wfs_subaper_offset(2) + wfs_subaper_space*(yc-1);
   end
end

if (dm_slave_corners_three_ways == 1)
   comp_nmsmax = 9;  % potentially we can have 8 masters, mastering some slave
else
   comp_nmsmax = 5;  % potentially we can have 4 masters, mastering some slave
end
comp_slave_defs = zeros(comp_nslav,comp_nmsmax);
comp_slave_weights = zeros(comp_nslav,comp_nmsmax);

% We will use this array to store the distances temporarily, so that
% we can calculate the weights properly.
comp_dist2 = zeros(comp_nmsmax - 1);

nearest1 = [1 0 -1 0; 0 1 0 -1];
nearest2 = [1 -1 -1 1; 1 1 -1 -1];
ic = 0;
nc = 0;
acount = 0;

for yc = 1:dm_act_num
   for xc = 1:dm_act_num
      if (comp_mx_acttype(yc,xc) == 2)
         nc = nc + 1;
         acount = 0;
         for ic = 1:4
            nxval = nearest1(1,ic) + xc;
            nyval = nearest1(2,ic) + yc;
            if ((nxval >= 1) & (nxval <= dm_act_num) & ...
                (nyval >= 1) & (nyval <= dm_act_num))
               if (comp_mx_acttype(nyval,nxval) == 1)
                  acount = acount + 1;
                  comp_slave_defs(nc,acount) = comp_mx_actnum(nyval,nxval);
               end
            end
         end
         near1 = acount;
         if ((acount == 0) | (dm_slave_corners_three_ways == 1))
            for ic = 1:4
               nxval = nearest2(1,ic) + xc;
               nyval = nearest2(2,ic) + yc;
               if ((nxval >= 1) & (nxval <= dm_act_num) & ...
                   (nyval >= 1) & (nyval <= dm_act_num))                   
                  if (comp_mx_acttype(nyval,nxval) == 1)
                     acount = acount + 1;
                     comp_slave_defs(nc,acount) = comp_mx_actnum(nyval,nxval);
                  end
               end
            end        
         end
         near2 = acount - near1;
         if ((dm_slave_corners_three_ways == 0) | (near1 == 0) | (near2 == 0))
            for ic = 1:acount          
               comp_slave_weights(nc,ic) = 1/acount;
            end
         else
            % At this point we know that we have near1 masters that are distance 1
            % away from current slave and near2 that are SQRT(2) away from it.
            if ((near1 > 0) & (near2 > 0))
               weight_near1 = sqrt(2)/(near1*sqrt(2) + near2);
               weight_near2 = 1/(near1*sqrt(2) + near2);
               for ic = 1:near1
                  comp_slave_weights(nc,ic) = weight_near1;
               end
               for ic = (near1+1):near2
                  comp_slave_weights(nc,ic) = weight_near2;
               end
            end
         end
      end
   end
end
