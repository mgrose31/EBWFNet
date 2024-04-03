function dmplot(doCompute)
% SYNTAX:
% dmplot(doCompute)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% doCompute [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: dmplot.m 3557 2015-10-25 00:16:14Z kbeardmore $

%% BEGIN_CODE

dmglobal;

reset_ax;  %reset both axis

if (doCompute == 1)
  dmcomp;    % here we compute different stuff.
end

d_subap;   % draw subapertures.  
plotcirc;  % first we draw all circles

slavnum = 0;
ic = 0;

h_actuators = zeros((dm_act_num * dm_act_num), 3);

h_slave2masterLines = zeros(comp_nslav, comp_nmsmax);
h_Lines = zeros(comp_nslav * comp_nmsmax, 4);
line_num = 1;

current_actuator_number = 1;
for yc = 1:dm_act_num
   for xc = 1:dm_act_num
   if (comp_mx_acttype(yc,xc) == 1)
      pntcol = 'white';
   elseif (comp_mx_acttype(yc,xc) == 2)
      pntcol = 'magenta';
   elseif (comp_mx_acttype(yc,xc) == 3)
      pntcol = 'green';
   end
      if (op_disp_slaves & (comp_mx_acttype(yc,xc) == 2))
         slavnum = slavnum + 1;
         for ic = 1:comp_nmsmax
            actvalx = comp_slave_defs(slavnum,ic);
            actvaly = comp_slave_defs(slavnum,ic);
            if ((actvalx > 0) & (actvaly > 0))
               actvalx = comp_vec_actnum(actvalx,1);
               actvaly = comp_vec_actnum(actvaly,2);
               xpoints = [comp_mx_actx(actvaly,actvalx), comp_mx_actx(yc,xc)];
               ypoints = [comp_mx_acty(actvaly,actvalx), comp_mx_acty(yc,xc)];
   
               tmp_Dx = abs(xpoints(2) - xpoints(1));
               tmp_Dy = abs(ypoints(2) - ypoints(1));
               tmp_D  = dm_act_spacing(1)*0.16;

               d_x = tmp_D * tmp_Dx / (tmp_Dx + tmp_Dy);
               d_y = tmp_D * tmp_Dy / (tmp_Dx + tmp_Dy);

               if (xpoints(2) > xpoints(1))
                  xpoints(1) = xpoints(1) + d_x;
                  xpoints(2) = xpoints(2) - d_x;
               else
                  xpoints(1) = xpoints(1) - d_x;
                  xpoints(2) = xpoints(2) + d_x;
               end
               if (ypoints(2) > ypoints(1))
                  ypoints(1) = ypoints(1) + d_y;
                  ypoints(2) = ypoints(2) - d_y;
               else
                  ypoints(1) = ypoints(1) - d_y;
                  ypoints(2) = ypoints(2) + d_y;
               end
               h_slave2masterLines(slavnum, ic+1) = patch('XData', xpoints,...
                                                          'YData', ypoints,...
                                                          'EdgeColor', 'yellow',...
                                                          'FaceColor', 'yellow');
               h_Lines(line_num, 1) = h_slave2masterLines(slavnum, ic+1);
               h_Lines(line_num, 2) = comp_mx_actnum(yc,xc);           % # of Slave act.
               h_Lines(line_num, 3) = comp_mx_actnum(actvaly,actvalx); % # of Master act.
               h_Lines(line_num, 4) = slavnum;
               line_num = line_num + 1;
            end
         end
      end
      if (op_disp_nums & (comp_mx_acttype(yc,xc) > 0))
         h_actuators(current_actuator_number, 2) = xc;
         h_actuators(current_actuator_number, 3) = yc;
         xpoints = comp_mx_actx(yc,xc);
         ypoints = comp_mx_acty(yc,xc);
         h_actuators(current_actuator_number, 1) = text('position',...
                                                        [xpoints,ypoints],'color',pntcol,...
                                                        'string',num2str(comp_mx_actnum(yc,xc)),...
                                                        'horizontalalignment','center',...
                                                        'verticalalignment','middle', 'EraseMode', 'xor');
      elseif (comp_mx_acttype(yc,xc) > 0)
         h_actuators(current_actuator_number, 2) = xc;
         h_actuators(current_actuator_number, 3) = yc;
         xpoints = [1 0 -1 0];
         ypoints = [0 1 0 -1];
         xpoints = xpoints*dm_act_spacing(1)*0.16 + comp_mx_actx(yc,xc);
         ypoints = ypoints*dm_act_spacing(2)*0.16 + comp_mx_acty(yc,xc);
         h_actuators(current_actuator_number, 1) = patch('XData', xpoints, ...
                                                         'YData', ypoints,...
                                                         'EdgeColor', pntcol,... 
                                                         'FaceColor', pntcol);
                                                         %'EraseMode', 'normal');
                                                         %'EraseMode', 'xor');
      end
      if ((comp_mx_acttype(yc,xc) == 2) & op_disp_slaves)
         h_slave2masterLines(slavnum, 1) = h_actuators(current_actuator_number, 1);
      end
      current_actuator_number = current_actuator_number + 1;
   end
end

% Now we draw the counts, if necessary:

if (op_disp_counts)   
   legend_x = dm_act_radii(3,2) + 0.05;
   legend_y = dm_act_radii(3,2) - 0.05;
   %draw # of inert act.
   dist_x = dm_act_spacing(1)*0.16*2;
   ax = gca;
   yl = get(gca,'ylim');
   yl = yl(2) - yl(1);
   dist_y = yl / 15; % dm_act_spacing(2)*0.16*3;
   comp_ninert = comp_nact - (comp_nmas + comp_nslav);
   xpoints = [1 0 -1 0];
   ypoints = [0 1 0 -1];
   xpoints = xpoints*dm_act_spacing(1)*0.16 + legend_x;
   ypoints = ypoints*dm_act_spacing(2)*0.16 + legend_y;
   text_xpoints = legend_x + 2*dist_x;
   text_ypoints = legend_y;
   patch(xpoints,ypoints,'green');
   text('position',[text_xpoints,text_ypoints],...
        'color','green',...
        'string',[num2str(comp_ninert) ' inert act'],...
        'horizontalalignment','left',...
        'verticalalignment','middle');
   ypoints = ypoints - dist_y;
   text_ypoints = text_ypoints - dist_y;
   patch(xpoints,ypoints,'magenta');
   text('position',[text_xpoints,text_ypoints],...
        'color','yellow',...
        'string',[num2str(comp_nslav) ' slave act'],...
        'horizontalalignment','left',...
        'verticalalignment','middle');
   ypoints = ypoints - dist_y;
   text_ypoints = text_ypoints - dist_y;
   patch(xpoints,ypoints,'white');
   text('position',[text_xpoints,text_ypoints],...
        'color','white',...
        'string',[num2str(comp_nmas) ' master act'],...
        'horizontalalignment','left',...
        'verticalalignment','middle');
   xpoints = [-1 1  1 -1 -1];
   ypoints = [ 1 1 -1 -1  1];
   xpoints = xpoints*dm_act_spacing(1)*0.16 + legend_x;
   ypoints = ypoints*dm_act_spacing(2)*0.16 + legend_y - dist_y*3;
   text_ypoints = text_ypoints - dist_y;
   line(xpoints, ypoints, 'color', 'red');
   % patch(xpoints,ypoints,'red');
   text('position',[text_xpoints,text_ypoints],...
        'color','red',...
        'string',[num2str(comp_nsubap) ' subap'],...
        'horizontalalignment','left',...
        'verticalalignment','middle');
end
