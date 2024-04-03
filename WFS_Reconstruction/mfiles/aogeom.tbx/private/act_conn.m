function act_conn(slaveActHandle, masterActHandle)
% SYNTAX:
% act_conn(slaveActHandle, masterActHandle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: act_conn.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE

dmglobal;
     
 res = find(h_actuators(:, 1) == masterActHandle);
 if(isempty(res) == 0)
   act_x_ind_end = h_actuators(res(1), 2);
   act_y_ind_end = h_actuators(res(1), 3);
 else
   return; % currentObject is not an actuator.
 end   

 res1 = find(h_actuators(:, 1) == slaveActHandle);
 if(isempty(res1) == 0)
   act_x_ind_start = h_actuators(res1(1), 2);
   act_y_ind_start = h_actuators(res1(1), 3);
 end
   
%fprintf(1, 'Trying to connect slave(%d) to master(%d) ...\n',...
%	comp_mx_actnum(act_y_ind_start, act_x_ind_start),...
%	comp_mx_actnum(act_y_ind_end, act_x_ind_end));



res2 = find(h_slave2masterLines(:, 1) == slaveActHandle);
if(isempty(res2) == 0)
  slav_num = res2(1);
end

res21 = find(comp_slave_defs(slav_num, 1:(comp_nmsmax-1)) == ...
			comp_mx_actnum(act_y_ind_end, act_x_ind_end));

if(isempty(res21) == 0) % Connection already exists
   return;		% So we do nothing,  just return.
end


res3 = find(comp_slave_defs(slav_num, 1:(comp_nmsmax-1)) == 0);

if(isempty(res3) == 1) % We need to add a column to a matrix.

  comp_nmsmax = comp_nmsmax + 1;
  %Adding rows to all 3 matrices:
  s = size(comp_slave_defs); % [#Rows, #NCols]
  comp_slave_defs    = [comp_slave_defs'; zeros(1, s(1))]';
  comp_slave_weights = [comp_slave_weights'; zeros(1, s(1))]';
  h_slave2masterLines= [h_slave2masterLines'; zeros(1, s(1))]';

  res3 = [(comp_nmsmax -1)];
end

  % creating a line and Recording new connection info:

  comp_slave_defs(slav_num, res3(1)) = comp_mx_actnum(act_y_ind_end, act_x_ind_end);

  actvalx = comp_slave_defs(slav_num, res3(1));
  actvaly = comp_slave_defs(slav_num, res3(1));
  
  if ((actvalx > 0) & (actvaly > 0))
    actvalx = comp_vec_actnum(actvalx,1);
    actvaly = comp_vec_actnum(actvaly,2);
    slave_act_num = comp_mx_actnum(act_y_ind_start, act_x_ind_start);
    xc = comp_vec_actnum(slave_act_num, 1);
    yc = comp_vec_actnum(slave_act_num, 2);

    xpoints = [comp_mx_actx(actvaly,actvalx), comp_mx_actx(yc,xc)];
    ypoints = [comp_mx_acty(actvaly,actvalx), comp_mx_acty(yc,xc)];

    tmp_Dx = abs(xpoints(2) - xpoints(1));
    tmp_Dy = abs(ypoints(2) - ypoints(1));
    tmp_D  = dm_act_spacing(1)*0.16;

    d_x = tmp_D * tmp_Dx / (tmp_Dx + tmp_Dy);
    d_y = tmp_D * tmp_Dy / (tmp_Dx + tmp_Dy);

    if(xpoints(2) > xpoints(1))
       xpoints(1) = xpoints(1) + d_x;
       xpoints(2) = xpoints(2) - d_x;
    else
       xpoints(1) = xpoints(1) - d_x;
       xpoints(2) = xpoints(2) + d_x;
    end
    if(ypoints(2) > ypoints(1))
       ypoints(1) = ypoints(1) + d_y;
       ypoints(2) = ypoints(2) - d_y;
    else
       ypoints(1) = ypoints(1) - d_y;
       ypoints(2) = ypoints(2) + d_y;
    end
  end

  h_slave2masterLines(slav_num, (res3(1)+1)) = patch('XData', xpoints,...
						   'YData', ypoints,...
						   'EdgeColor', 'yellow',...
						   'FaceColor', 'yellow');

  % Adding a row to a h_Lines:

  res4 = find(h_Lines(:, 1) == 0);
  if(isempty(res4) == 0)
    line_num = res4(1);
  else
    h_Lines  = [h_Lines; zeros(1, 4)];
    s        = size(h_Lines);
    line_num = s(1);
  end

  h_Lines(line_num, 1) = h_slave2masterLines(slav_num, (res3(1) + 1));
  h_Lines(line_num, 2) = comp_mx_actnum(act_y_ind_start, act_x_ind_start); % Slave
  h_Lines(line_num, 3) = comp_mx_actnum(act_y_ind_end, act_x_ind_end);     % Master
  h_Lines(line_num, 4) = slav_num;
