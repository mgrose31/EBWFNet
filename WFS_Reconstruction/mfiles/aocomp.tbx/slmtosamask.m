function samask = slmtosamask(wg, slmmask)
% SYNTAX: 
% samask = slmtosamask(wg, slmmask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: slmtosamask.m 3027 2010-09-21 21:04:10Z amoran $

%% BEGIN_CODE

samask = zeros(wg.nx, wg.ny);
for jj=1:wg.ny
   for ii=1:wg.nx
      if (slmmask(round(wg.x(ii)),round(wg.y(jj))))
         samask(ii,jj) = 1;
      end;
   end;
end;
