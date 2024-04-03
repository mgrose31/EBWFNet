function connect(currentObject)
% SYNTAX:
% connect(currentObject)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: connect.m 3028 2010-09-21 21:04:58Z amoran $

%% BEGIN_CODE

dmglobal;

%At this point we have lastSelectedObject and currentObject

 res = find(h_actuators(:, 1) == currentObject);
 if(isempty(res) == 0)
   act_x_ind_end = h_actuators(res(1), 2);
   act_y_ind_end = h_actuators(res(1), 3);
 else
   return; % currentObject is not an actuator.
 end   

 res1 = find(h_actuators(:, 1) == lastSelectedObject(1));
 if(isempty(res1) == 0)
   act_x_ind_start = h_actuators(res1(1), 2);
   act_y_ind_start = h_actuators(res1(1), 3);
 end
   
 if((comp_mx_acttype(act_y_ind_start, act_x_ind_start) == 2) &...% Slave
    (comp_mx_acttype(act_y_ind_end, act_x_ind_end) == 1))	% Master
   act_conn(lastSelectedObject(1), currentObject);  % act_conn(Slave, Master)
 elseif((comp_mx_acttype(act_y_ind_start, act_x_ind_start) == 1) &... % Master
    (comp_mx_acttype(act_y_ind_end, act_x_ind_end) == 2))	   % Slave
   act_conn(currentObject, lastSelectedObject(1));  % act_conn(Slave, Master)
 end
   