function [imin, imax, simg, ncol, imghndl] = tvscl(img, x, y)
% SYNTAX:
% [imin, imax, simg, ncol, imghndl] = tvscl(img, x, y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% img [ ] = 
% x [ ] =
% y [ ] =
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% imin [ ] =
% imax [ ] =
% simg [ ] =
% ncol [ ] =
% imghndl [handle] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: tvscl.m 3051 2010-10-01 20:33:26Z amoran $

%% BEGIN_CODE


si = size(img);
if ((si(1) == 1) | (si(2) == 1))
   n = si(1)*si(2);
   nr = round(sqrt(n));
   if ((nr*nr) == n)
      img = reshape(img, nr, nr);
   end;
end;
imin = min(min(img));
imax = max(max(img));
map = colormap;
sm = size(map);
ncol = sm(1);
simg = ((img - imin)./(imax - imin)) .* ncol;
if (exist('x') & exist('y'))
   imghndl = image(x, y, simg);
else
   imghndl = image(simg);
end;
