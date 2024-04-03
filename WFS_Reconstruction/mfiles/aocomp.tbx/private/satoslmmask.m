function slmmask = satoslmmask(wg, samask)
% SYNTAX: 
% slmmask = satoslmmask(wg, samask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: satoslmmask.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

slmmask = zeros(wg.xmax, wg.ymax);
for jj=1:wg.ny
   for ii=1:wg.nx
      if (samask(ii,jj))
         slmmask(wg.u(ii):wg.d(ii),wg.l(jj):wg.r(jj)) = 1;
      end;
   end;
end;
