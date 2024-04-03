function unselect(object)
% SYNTAX:
% unselect(object)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: unselect.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE

if(object(4) > 0) 
  delete(object(4));
  object(4) = -1;
  object(1) = -1;
end
