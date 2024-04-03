function dmbdown
% SYNTAX:
% dmbdown
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: dmbdown.m 3028 2010-09-21 21:04:58Z amoran $

%% BEGIN_CODE

dmglobal;

isMouseButtonDown = 1;

set(gui_main, 'Units', 'normalized');
currentPoint  = get(gui_main, 'CurrentPoint');
set(gui_main, 'Units', 'pixels');
selectionType = get(gui_main, 'SelectionType');

currentObject = get(gui_main, 'CurrentObject');
currentTime   = clock;

if((lastSelectedObject(2) == ACTUATOR) & (strcmp(selectionType, 'alt')==1) &...
   (op_disp_slaves == 1))
     connect(currentObject);
     unselect(lastSelectedObject);
     lastClickTime = -1;
     lastSelectedObject = [-1 -1 -1 -1];
     return;
end

if(lastClickTime == -1)
  lastSelectedObject = sel_obj(currentObject);
  lastClickTime = currentTime;
else
%lastSelectedObject
%currentObject
% (currentTime(6) - lastClickTime(6))
  if(((lastSelectedObject(1)==currentObject) | (lastSelectedObject(4)==currentObject)) & ...
     (lastSelectedObject(3) == SUPPORTS_DOUBLE_CLICK) & ...
     (currentTime(6) - lastClickTime(6) < doubleClickInterval))
    activate(lastSelectedObject);
    unselect(lastSelectedObject);
    lastSelectedObject = [-1 -1 -1 -1];
  else
    if((lastSelectedObject(1) == currentObject) | (lastSelectedObject(4) == currentObject))
      unselect(lastSelectedObject);
      lastSelectedObject = [-1 -1 -1 -1];
    else
      if(lastSelectedObject(1) > 0) 
        unselect(lastSelectedObject);
      end;
      lastSelectedObject = sel_obj(currentObject);
    end
    lastClickTime = currentTime;
  end
end
