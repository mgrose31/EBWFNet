function dmbup
% SYNTAX:
% dmbup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: dmbup.m 3028 2010-09-21 21:04:58Z amoran $

%% BEGIN_CODE

dmglobal;

isMouseButtonDown = 0;
objUp = get(gui_main, 'CurrentObject');

