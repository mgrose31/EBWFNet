function z = circ(x, y, d)
% SYNTAX:
% z = circ(ny, nx, r)
%   where
%      ny = number of points in y
%      nx = number of points in x
%      r = radius in number of points
%           OR
% z = circ(x,y,d)
%   where
%      x = mesh of x coordinates
%      y = mesh of y coordinates
%      d = diameter in coordinate space
%   x and y are created via meshgrid.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: circ.m 4212 2023-05-05 17:52:03Z jtellez $

%% BEGIN_CODE

narginchk(2,3);
if nargin<3, d=1; end

if isscalar(x) && isscalar(y)
    z = ones(y, x);
    r = d;
    ny2 = y/2;
    nx2 = x/2;
    rsq = r.^2;
    for k=1:x
        ksq = (k-nx2).^2;
        if (ksq <= rsq)
            for j=1:y
                if ((((j-ny2).^2)+ksq) >= rsq)
                    z(j,k) = 0;
                end
            end
        else
            z(:,k) = 0;
        end
    end
else
    r = sqrt(x.*x+y.*y);
    z = double(r<d/2);
    z(r==d/2)=0.5;
end