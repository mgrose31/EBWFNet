function dmwindon(windowname)
% SYNTAX:
% dmwindon(windowname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: dmwindon.m 3359 2011-02-25 15:58:12Z morris $

%% BEGIN_CODE

dmglobal;

if (strcmp(windowname,'gui_rectangle'))
   set(gui_rectangle,'visible','on');
else
   eval(windowname);
end
