function result = sel_obj(currentObject)
% SYNTAX:
% result = sel_obj(currentObject)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: sel_obj.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE

dmglobal;

     
  % First we check if our current object is an actuator
  % because an actuator is the smallest object and it is in the
  % circle's hot area. This is stupid,  but Matlab doesn't let us
  % change the hot area of an object.   ... unles we write
  % our own 'inside' function,  which of course will be much slower.

% 1st element - handle of the selected object
% 2nd element - type of an object
% 3rd element - whether it supports doubleclicking
% 4th element - handle of the selectedBorder, or -1 if no border is drawn
%

result = [-1 -1 -1 -1];

result(1) = currentObject;
result(3) =  ~SUPPORTS_DOUBLE_CLICK; % does not support doubleClicking;
result(4) = -1;                      % object that denotes selection, if any

res = find(h_actuators(:,1) == currentObject);
if(isempty(res) == 0)
%    fprintf('User clicked on an actuator at (%d, %d)\n',...
%             h_actuators(res(1), 3), h_actuators(res(1), 2));
    result(2) = ACTUATOR;    
    result(3) = SUPPORTS_DOUBLE_CLICK;    
    tmp_x = comp_mx_actx(h_actuators(res(1), 3), h_actuators(res(1), 2));
    tmp_y = comp_mx_acty(h_actuators(res(1), 3), h_actuators(res(1), 2));

    tmp_radius = max(dm_act_spacing)*0.18;
    result(4) = newcircl(tmp_radius, tmp_x, tmp_y, 'yellow');
    return;
end

res = find(h_Lines(:,1) == currentObject);
if(isempty(res) == 0)
%	fprintf(1, 'Selecting a S2M Line\n');
    result(2) = SLAVE_TO_MASTER_LINE;
    result(3) = SUPPORTS_DOUBLE_CLICK;
    xpts = get(result(1), 'XData');
    ypts = get(result(1), 'YData');
    result(4) = line(xpts, ypts, 'color', 'green', 'EraseMode', 'xor');
    return;
end


if(currentObject == h_inertCircles(1)) 
%   fprintf('User clicked on a inner Inert Circle\n');
    result(2) = INERT_INNER_CIRCLE;    
elseif(currentObject == h_inertCircles(2)) 
%   fprintf('User clicked on a outer Inert Circle\n');
    result(2) = INERT_OUTER_CIRCLE;    

elseif(currentObject == h_slaveCircles(1)) 
%   fprintf('User clicked on a inner Slave Circle\n');
   result(2) = SLAVE_INNER_CIRCLE;    
elseif(currentObject == h_slaveCircles(2)) 
%   fprintf('User clicked on a outer Slave Circle\n');
   result(2) = SLAVE_OUTER_CIRCLE;    

elseif(currentObject == h_masterCircles(1)) 
%   fprintf('User clicked on a inner Master Circle\n');
   result(2) = MASTER_INNER_CIRCLE;    
elseif(currentObject == h_masterCircles(2)) 
%   fprintf('User clicked on a outer Master Circle\n');
   result(2) = MASTER_OUTER_CIRCLE;    

elseif(currentObject == h_clampCircle) 
%   fprintf('User clicked on a clamp Circle\n');
   result(2) = CLAMPING_CIRCLE;
elseif(currentObject == h_subapClipCircle) 
%   fprintf('User clicked on a subapClipCircle\n');
   result(2) = SUBAPERTURE_OUTER_CIRCLE;
elseif(currentObject == h_subapAnulusCircle) 
%   fprintf('User clicked on a subaperture Anulus Circle\n');
   result(2) = SUBAPERTURE_INNER_CIRCLE;
else
%   fprintf('User clicked on some unknown object\n');
   result(1) = -1; 
end

