function dmsave(filename)
% SYNTAX:
% dmsave(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% filename [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% handle [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: dmsave.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE

dmglobal;

nact = comp_nact;
nslav = comp_nslav;
nmas = comp_nmas;

xact = zeros(nact,1);
yact = zeros(nact,1);
actype = zeros(nact,1);
u = zeros(nact,1);
w = zeros(nact,1);

mx_u = zeros(dm_act_num,dm_act_num);
mx_w = zeros(dm_act_num,dm_act_num);
centerval = floor(dm_act_num/2) + 1;
mx_u(centerval,centerval) = 200;
mx_w((centerval-1):(centerval+1),(centerval-1):(centerval+1)) = ...
     [0.01 0.05 0.01; 0.05 1 0.05; 0.01 0.05 0.01] .* 3.3E-6;
   
ic = 0;
for yc = 1:dm_act_num
   for xc = 1:dm_act_num
      if (comp_mx_acttype(yc,xc) > 0)
         ic = ic + 1;
         xact(ic) = comp_mx_actx(yc,xc);
         yact(ic) = comp_mx_acty(yc,xc);
         if (comp_mx_acttype(yc,xc) == 1)
            actype(ic) = 1;
         elseif (comp_mx_acttype(yc,xc) == 2)
            actype(ic) = -1;
         else
            actype(ic) = 0;
         end
         u(ic) = mx_u(yc,xc);
         w(ic) = mx_w(yc,xc);
      end
   end
end

%This is what it is now:
nmsmax = 1;
for i = 1:(comp_nmsmax - 1)
   if (any(comp_slave_defs(:, i)) == 1)
      nmsmax = nmsmax + 1;
   end
end

% At this point we have the number of columns that we need to store in
% the file.
idxmas = zeros(comp_nslav, nmsmax);
weight = zeros(comp_nslav, nmsmax);
cnt = 1;
for i = 1:(nmsmax-1)
   if (any(comp_slave_defs(:, i)) == 1)
      idxmas(:,cnt) = comp_slave_defs(:, i);
      weight(:,cnt) = comp_slave_weights(:, i);
      cnt = cnt + 1;
   end
end % comp_slave_defs

%Now we need to compute weights for each row :
for i = 1:comp_nslav
   res = find(h_Lines(:, 4) == i);
   if (isempty(res) == 0)
      s = size(res);
      d = zeros(1, s(1));
      lines = h_Lines(res, 1);
      for j = 1:s(1)
         line = lines(j);
         xpts = get(line, 'XData');
         ypts = get(line, 'YData');
         % Distance between slave and master    
         d(j) = sqrt((xpts(2) - xpts(1))^2 + (ypts(2) - ypts(1))^2);
      end
      ww = (1 ./ d) / sum(1 ./ d); % We construct a vector of weights from
      k = 1; % the vector of distances.
      for j = 1:(nmsmax - 1)
         if (idxmas(i, j) > 0)
            weight(i, j) = ww(k);
            k = k + 1;
         else
            weight(i, j) = 0;
         end
      end
   end
end

xsub = 0;
ysub = 0;
nsub = 0;

subaprad = 0.5*wfs_full_only*sqrt(wfs_subaper_space^2+wfs_subaper_space^2);
for xc = 1:wfs_act_num
   for yc = 1:wfs_act_num
      d = ((comp_mx_wfsx(yc,xc) - wfs_subaper_offset(1))^2 + ...
           (comp_mx_wfsy(yc,xc) - wfs_subaper_offset(2))^2);
      if ((d <= (wfs_subaper_clip - subaprad)^2) & ...
          (d >= (wfs_subaper_anulus + subaprad)^2))
         nsub = nsub + 1;
         xsub(nsub) = comp_mx_wfsx(yc,xc);
         ysub(nsub) = comp_mx_wfsy(yc,xc);
      end
   end
end

xsub = xsub';
ysub = ysub';

actnum = dm_act_num;
subnum = wfs_act_num;
rc = dm_clamp_rad;
ninfl = dm_ninfl;
nsi = wfs_num_steps;
hs = wfs_subaper_space;
rmas_in  = dm_act_radii(1,1);
rslav_in = dm_act_radii(2,1);
ract_in  = dm_act_radii(3,1);
rsub_in  = wfs_subaper_anulus;

rmas_out  = dm_act_radii(1,2);
rslav_out = dm_act_radii(2,2);
ract_out  = dm_act_radii(3,2);
rsub_out  = wfs_subaper_clip;

xoff = dm_act_offset(1);
yoff = dm_act_offset(2);
center_xoff = dm_act_center_offset(1);
center_yoff = dm_act_center_offset(2);
slave_corners_3ways = dm_slave_corners_three_ways;

xoffsub = wfs_subaper_offset(1);
yoffsub = wfs_subaper_offset(2);
fullsub = wfs_full_only;

dminfl = 0;
if (dm_influence == 1); dminfl = 0; end;
if (dm_influence == 2); dminfl = 4; end;
if (dm_influence == 3); dminfl = 1; end;
if (dm_influence == 4); dminfl = 2; end;
if (dm_influence == 5); dminfl = 3; end;
xwidth = dm_infl_width(1);
ywidth = dm_infl_width(2);

% MLM add 16 Jan 2007 - try to match expecation of TasatDMModel in
% AOIInf.cpp. I am basing the value saved on a comment in TasatDMModel.cpp
% which states that "For Green's inf. fun. Bal's model is used.". balmodel 
% appears to be a flag, so if the selected influence function is Green's, 
% I set balmodel to 1, else to 0;
%
% RWPII mod 20080530. The comment in TasatDMModel.cpp is ambiguous and
% could reasonably lead someone to believe that the above comment is the
% right course of action. But it is not. The balmodel=1 option is vestigal
% and should only be used by someone who really knows that he wants it.
% Moreover, it is such a special case, that always setting balmodel=0
% here is the right thing to do. If a user wants to use balmodel=1, he can
% make a hand mod. The Greens function influence function is not dependent
% on balmodel=1. Rather the balmodel=1 option handles an actuator cross-
% talk issue that exists for certain deformable mirror types which are
% not presently used a lot.
%
% So in a nutshell, always set balmodel=0.
%
balmodel = 0;
%if (dminfl == 0) balmodel = 1; end; 

save(filename,'nact','nslav','nmas','xact','yact','actype',...
     'nmsmax','idxmas','weight','nsub','xsub','ysub','u','w',...
     'rc','ninfl','nsi','hs','rmas_in','rslav_in','ract_in','rsub_in',...
     'rmas_out','rslav_out','ract_out','rsub_out',...
     'actnum','subnum','xoff','yoff','center_xoff', 'center_yoff',...
     'slave_corners_3ways', 'xoffsub','yoffsub',...
     'fullsub','dminfl','xwidth','ywidth','balmodel');