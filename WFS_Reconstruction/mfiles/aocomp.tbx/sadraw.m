function sadraw(wg, marks, mask)
% SYNTAX:
% sadraw(wg, marks, mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% wg [ ] = 
% marks [ ] = 
% mask [ ] =
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% sadraw [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: sadraw.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

nx = prod(size(wg.x));
ny = prod(size(wg.y));
if (~exist('mask','var'))
   mask = ones(nx,ny);
end;
if (isempty(mask))
   mask = ones(nx,ny);
end;
if (~ishold), hold; end;
for ii=1:ny
   for jj=1:nx
      if (mask(jj,ii))
         plot([wg.l(ii),wg.r(ii),wg.r(ii),wg.l(ii),wg.l(ii)],...
              [wg.u(jj),wg.u(jj),wg.d(jj),wg.d(jj),wg.u(jj)],'w');
      end;
   end;
end;
if (~exist('marks','var'))
   marks = [];
end;
if (~isempty(marks))
   for ii=1:ny
      for jj=1:nx
         if (mask(jj,ii))
            plot(wg.y(ii)+marks(1,ii,jj), wg.x(jj)+marks(2,ii,jj), 'r+');
         end;
      end;
   end;
end;
hold;