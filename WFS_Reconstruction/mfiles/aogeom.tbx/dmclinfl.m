function dmclinfl
% SYNTAX:
% dmclinfl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: dmclinfl.m 3028 2010-09-21 21:04:58Z amoran $

%% BEGIN_CODE

dmglobal;

myobj = get(gcf,'currentobject');
temp_dminf = 0;
for ic = 1:5
   if (guii_function(ic) == myobj)
      set(myobj,'value',1);
      temp_dminf = ic;
   else
      set(guii_function(ic),'value',0);
   end
end

tempstr = get(myobj,'userdata');
set(guii_infl_label(1),'string', tempstr(1,:));
set(guii_infl_label(2),'string', tempstr(2,:));

if (temp_dminf == 1)
   set(guii_infl_width(1),'string',num2str(dm_clamp_rad)); 
   set(guii_infl_width(2),'string',num2str(dm_ninfl));
else
   set(guii_infl_width(1),'string',num2str(dm_infl_width(1))); 
   set(guii_infl_width(2),'string',num2str(dm_infl_width(2)));
end
