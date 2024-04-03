function [opd] = dmmap(ixact, iyact, xopd, yopd, opdif, dact);
% SYNTAX:
% [opd] = dmmap(ixact, iyact, xopd, yopd, opdif, dact);
% [opd] = dmmap(wsdm, xopd, yopd, dact, method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% ixact [ ] = 
% iyact [ ] = 
% xopd [ ] = 
% yopd [ ] = 
% opdif [ ] = 
% dact [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% opd [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: dmmap.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

if (nargin == 6)
   nx = prod(size(xopd));
   ny = prod(size(yopd));
   opd = zeros(nx,ny);
   nact = prod(size(ixact));
   so = size(opdif);
   nside = sqrt(so(2));
   nhalf = floor(nside/2);
   for i=1:nact
      subif = reshape(opdif(i,:), nside, nside) * dact(i);
      ix = ixact(i);
      jy = iyact(i);
      ilow = max([1, ix-nhalf]);
      ihigh = min([nx, ix+nhalf]);
      jlow = max([1, jy-nhalf]);
      jhigh = min([ny, jy+nhalf]);
      ilows = ilow - ix + nhalf + 1;
      ihighs = ihigh - ix + nhalf + 1;
      jlows = jlow - jy + nhalf + 1;
      jhighs = jhigh - jy + nhalf + 1;
      opd(ilow:ihigh,jlow:jhigh) = opd(ilow:ihigh,jlow:jhigh) + subif(ilows:ihighs,jlows:jhighs);
   end;
else
   dact = yopd;
   yopd = xopd;
   xopd = iyact;
   wsdm = ixact;
   if (nargin <= 4)
      method = '*linear';
   else
      method = opdif;
   end;
   if (~exist('method','var'))
      method = '*linear';
   end;
   if (isempty(method))
      method = '*linear';
   end;
   opd = zeros([length(xopd), length(yopd)]);
   hdxopd = 0.5*(xopd(2)-xopd(1));
   hdyopd = 0.5*(yopd(2)-yopd(1));
   for ii=1:length(wsdm.xact)
      if (iscell(wsdm.actinfl))
         ai.x = wsdm.xact(ii) + wsdm.actinfl{ii}.x;
         ai.y = wsdm.yact(ii) + wsdm.actinfl{ii}.y;
         ai.g = dact(ii) .* wsdm.actinfl{ii}.g;
      else
         ai.x = wsdm.xact(ii) + wsdm.actinfl.x;
         ai.y = wsdm.yact(ii) + wsdm.actinfl.y;
         ai.g = dact(ii) .* wsdm.actinfl.g;
      end;
      iind = find(((ai.x(1)-hdxopd) <= xopd) & (xopd <= (ai.x(end)+hdxopd)));
      if (~isempty(iind))
         jind = find(((ai.y(1)-hdyopd) <= yopd) & (yopd <= (ai.y(end)+hdyopd)));
         if (~isempty(jind))
            [aixx, aiyy] = meshgrid(ai.x, ai.y);
            [opdxx, opdyy] = meshgrid(xopd(iind), yopd(jind));
            dopd = interp2(aixx, aiyy, ai.g, opdxx, opdyy, method, 0)';
            opd(iind(1):iind(end),jind(1):jind(end)) = opd(iind(1):iind(end),jind(1):jind(end)) + dopd;
         end;            
      end;
   end;
end;  