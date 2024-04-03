function activate(object)
% SYNTAX:
% activate(object)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% object [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% handle [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: activate.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE

dmglobal;

if(object(2) == ACTUATOR)  

 res = find(h_actuators(:, 1) == object(1));
 if(isempty(res) == 0)
    %first we get indices of the selected actuator.
    act_x_ind = h_actuators(res(1), 2);
    act_y_ind = h_actuators(res(1), 3);
     
    if(comp_mx_acttype(act_y_ind, act_x_ind) == 1) % Master
      comp_mx_acttype(act_y_ind, act_x_ind) = 2;   % Set it to Slave
      if(op_disp_nums)
	set(h_actuators(res(1), 1), 'color', SlaveColor);
      else
	set(h_actuators(res(1), 1), 'EdgeColor', SlaveColor, 'FaceColor', SlaveColor);
      end

      act_num = comp_mx_actnum(act_y_ind, act_x_ind);

	% res2 is a vector that contains indices of the act we are looking for
	% as if the matrix comp_slave_defs were Nx1 matrix (columnwise)
      	res2 = find(comp_slave_defs == act_num);
	% So, to detect real indices we must do the following:
	csd_dim = size(comp_slave_defs);
	row_nums = rem(res2, csd_dim(1));
	col_nums = fix(res2 ./ csd_dim(1)) + 1;
	%At this point we have the coords of our master Actuator in comp_slave_defs matrix
	%and now we set those values to 0 and delete appropriate lines from the picture.
	size_res2 = size(res2);
	for i = 1:size_res2(1)
	  if(row_nums(i) == 0)
	    row_nums(i) = csd_dim(1);
	    col_nums(i) = col_nums(i) - 1;
	  end
	  comp_slave_defs(row_nums(i), col_nums(i)) = 0;
	  res3 = find(h_Lines(:, 1) == h_slave2masterLines(row_nums(i), (col_nums(i)+1)));
	  if(isempty(res3) == 0)
	    h_Lines(res3(1), :) = [];
	  end;
	  delete(h_slave2masterLines(row_nums(i), (col_nums(i)+1)));
	  h_slave2masterLines(row_nums(i), (col_nums(i)+1)) = 0;
	end
	% Now we need to add a row to both comp_slave_defs and comp_slave_weights matrices
	% as well as increment comp_nslav.
	num_of_slaves_before = 0;
	for yc = 1:dm_act_num
	   for xc = 1:dm_act_num
	     if((yc == act_y_ind) & (xc == act_x_ind))
		break;
	     elseif(comp_mx_acttype(yc, xc) == 2)
		num_of_slaves_before = num_of_slaves_before + 1;
	     end
	   end
     	   if((yc == act_y_ind) & (xc == act_x_ind)), break, end
	end
	% Now we need to insert rows into both matrices:
	s = size(comp_slave_defs);
	part1 = comp_slave_defs(1:num_of_slaves_before, :);
	part2 = comp_slave_defs((num_of_slaves_before+1):comp_nslav, :);
	comp_slave_defs = [part1; zeros(1, s(2)); part2]; % We insert row of zeros since
							  % our new Slave isn't connected
	part1 = comp_slave_weights(1:num_of_slaves_before, :);
	part2 = comp_slave_weights((num_of_slaves_before+1):comp_nslav, :);
	comp_slave_weights = [part1; zeros(1, s(2)); part2]; % We insert row of zeros since
							     % our new Slave isn't connected
	part1 = h_slave2masterLines(1:num_of_slaves_before, :);
	part2 = h_slave2masterLines((num_of_slaves_before+1):comp_nslav, :);
	h_slave2masterLines = [part1; zeros(1, s(2)); part2]; % We insert row of zeros since
							     % our new Slave isn't connected
	h_slave2masterLines((num_of_slaves_before+1), 1) = h_actuators(res(1), 1);
	comp_nslav = comp_nslav + 1;
    elseif(comp_mx_acttype(act_y_ind, act_x_ind) == 2) % Slave
      comp_mx_acttype(act_y_ind, act_x_ind) = 3;   % Set it to Inert
      if(op_disp_nums)
	set(h_actuators(res(1), 1), 'color', InertColor);
      else
	set(h_actuators(res(1), 1), 'EdgeColor', InertColor, 'FaceColor', InertColor);
      end
      res1 = find(h_slave2masterLines(:, 1) == h_actuators(res(1), 1));
      if(isempty(res1) == 0)
	%First we delete lines
	for i = 2:(comp_nmsmax)
	  if(h_slave2masterLines(res1(1), i) > 0)
	    res3 = find(h_Lines(:, 1) == h_slave2masterLines(res1(1), i));
	    if(isempty(res3) == 0)
	      h_Lines(res3(1), :) = [];
	    end;
	    delete(h_slave2masterLines(res1(1), i));
	  end
	end

	%Then we delete row in h_slave2masterLines
	h_slave2masterLines(res1(1),:) = [];

	%Finally, we delete row in comp_slave_defs and comp_slave_weights
	% We also need to decrement comp_nslav
	comp_slave_defs(res1(1), :) = [];
	comp_slave_weights(res1(1), :) = [];
	comp_nslav = comp_nslav - 1;
      end
    elseif(comp_mx_acttype(act_y_ind, act_x_ind) == 3) % Inert
      comp_mx_acttype(act_y_ind, act_x_ind) = 1;   % Set it to Master
      if(op_disp_nums)
	set(h_actuators(res(1), 1), 'color', MasterColor);
      else
	set(h_actuators(res(1), 1), 'EdgeColor', MasterColor, 'FaceColor', MasterColor);
      end
   end
 end

h_Lines = zeros(comp_nslav * comp_nmsmax, 4);
line_num = 1;
slavnum = 0;
for yc = 1:dm_act_num
    for xc = 1:dm_act_num
        if (op_disp_slaves & (comp_mx_acttype(yc,xc) == 2))
           slavnum = slavnum + 1;
           for ic = 1:comp_nmsmax
               actvalx = comp_slave_defs(slavnum,ic);
               actvaly = comp_slave_defs(slavnum,ic);
               if ((actvalx > 0) & (actvaly > 0))
                  actvalx = comp_vec_actnum(actvalx,1);
                  actvaly = comp_vec_actnum(actvaly,2);
		  h_Lines(line_num, 1) = h_slave2masterLines(slavnum, ic+1);
		  h_Lines(line_num, 2) = comp_mx_actnum(yc,xc);           % # of Slave act.
		  h_Lines(line_num, 3) = comp_mx_actnum(actvaly,actvalx); % # of Master act.
		  h_Lines(line_num, 4) = slavnum;
		  line_num = line_num + 1;
%		  if(line_num > 4), return, end
                end
            end
        end
     end
end
%h_Lines

elseif(object(2) == SLAVE_TO_MASTER_LINE)

%object
%h_Lines

  res = find(h_Lines(:, 1) == object(1));
  if(isempty(res) == 0)
	line_Handle    = h_Lines(res(1), 1)
	slave_act_num  = h_Lines(res(1), 2)
	master_act_num = h_Lines(res(1), 3)
	slave_num      = h_Lines(res(1), 4)
	
	% delete row in the h_Lines matrix:
	h_Lines(res(1), :) = [];
	
	res1 = find(comp_slave_defs(slave_num, :) == master_act_num);
	if(isempty(res1) == 0)
	  comp_slave_defs(slave_num, res1(1)) = 0;
	  comp_slave_weights(slave_num, res1(1)) = 0;
	  h_slave2masterLines(slave_num, (res1(1)+1)) = 0; % ( +1 since 1st el-t of row - 
							   %  a handle to a slave actuator)
	end
	delete(line_Handle);
  end
end




