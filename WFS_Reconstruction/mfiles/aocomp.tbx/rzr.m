function [anz] = rzr(a, releps);
% SYNTAX:
% [anz] = rzr(a, {releps})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% (r)emove (z)ero (r)ows. Remove all zero rows from the matrix a. When 
% releps is not specified a row is zero is all of its elements are 
% precisely zero. If releps is specified, then a row is zero if the 
% absolute value of all elemens are less than releps*max(max(a)).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% a [ ] = 
% releps [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% anz [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: rzr.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

sa = size(a);
anz = 0*a;
nz = 0;
if (nargin == 1)
   for k=1:sa(1)
      ind = find(a(k,:));
      if (prod(size(ind)) ~= 0)
         nz = nz + 1;
         anz(nz,:) = a(k,:);
      end;
   end;
   anz = anz(1:nz,:);
else
   epschk = releps*max(max(a));
   for k=1:sa(1)
      ind = find(abs(a(k,:)) >= epschk);
      if (prod(size(ind)) ~= 0)
         nz = nz + 1;
         anz(nz,:) = a(k,:);
      end;
   end;
   anz = anz(1:nz,:);
end;
