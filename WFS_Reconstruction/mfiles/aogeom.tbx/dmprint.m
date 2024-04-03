function dmprint(arg)
% SYNTAX:
% dmprint(arg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% arg [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: dmprint.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE

dmglobal;

figure(gui_main);
if(arg == 0)
  set(gcf, 'InvertHardCopy', 'on');
else
  set(gcf, 'InvertHardCopy', 'off');
end
print;
