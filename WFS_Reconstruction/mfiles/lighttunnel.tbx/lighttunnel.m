function [img, x, y] = lighttunnel(imgp, cx, cy, p, tx, ty, bx, by, reduce)
% SYNTAX:
% [img, x, y] = lighttunnel(imgp, cx, cy, p, tx, ty, bx, by, reduce)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Adds specified tilt and blur effects to input image.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% imgp [matrix] = Input image (must be double)
% cx [matrix] = X position of image pixels, increasing down columns (arb)
% cy [matrix] = Y position of image pixels, increasing across rows (arb)
% p [scalar/matrix] = power multiplier of output image values
% tx [matrix] = X-tilts for image pixel locations (pix)
% ty [matrix] = Y-tilts for image pixel locations (pix)
% bx [scalar/matrix] = Row blur function per pixel (pix)
% by [scalar/matrix] = Column blur function per pixel (pix)
% reduce = Option to reduce the output image with valid inputs:
%             0 = no reduction
%             1 = reduce image size to input image size
%             2 = reduce border further
%             3 = reduce border and return image as uint8
%          For rectangular images, reduce is set to 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% img [matrix] = imgp with applicable tilt and blur added image
% x [vector] = Image x-coordinates (increasing down row)
% y [vector] = Image y-coordinates (increasing across column)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEE ALSO: ltbound.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Robert W. Praus, II
% (c) 2005 MZA Associates Corporation, Albuquerque, NM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: lighttunnel.m 3036 2010-09-23 21:22:50Z amoran $


%% BEGIN_CODE


%% Number of Pixels in Each Dimension
[nx,ny,nz] = size(imgp);

%% Set Variable 'reduce'
if (~exist('reduce','var'))
    reduce = [];
end;

if (isempty(reduce))
    reduce = 1;
end;

if nx~=ny;
    if(reduce~=0)
        disp('Input image not square.  Reduction option off.')
        reduce=0;
    end
end

%% Call to ltBound
[lx, rx, minx, maxx, cxx] = ltbound(cx, tx, bx);
[ly, ry, miny, maxy, cyy] = ltbound(cy, ty, by);
% lx (or ly) is tx (tilt function per pixel)  - 0.5* the blurring function
% rx (or ry) is tx (tilt function per pixel) + 0.5* the blurring function
% minx is image left boundary
% maxx is image right boundary

%% Create Distorted Image, zeros plus border in order to output entire
%% distorted image
imgNrow = nx+(maxx-minx);
imgNcol = ny+(maxy-miny);
img = zeros(imgNrow,imgNcol,nz);

%% Iterates by column
for jj=miny:maxy        % Iterates by column
    ply = max(cyy+jj-0.5*bx, ly);   % Select maximum:  spatial y-coordinates + iterand - 0.5 * blurring, or input distortion function
    pry = min(cyy+jj+0.5*by, ry);   % Select minimum:  spatial y-coordinates + iterand + 0.5 * blurring, or input distortion function
    pdy = max(pry-ply, 0);          % thresholds pry-ply value at zero
    raty = pdy ./ by;               % Ratio of thresholded tilt/blur to blur function for this column
    for ii=minx:maxx                % Iterates by row (one iteration for each pixel in column)
        plx = max(cxx+ii-0.5*bx, lx);
        prx = min(cxx+ii+0.5*by, rx);
        pdx = max(prx-plx,0);
        ratx = pdx ./ bx;
        area = ratx .* raty;        % Thresholded tilt/blur to blur ratio in x (row) .* thresholded tilt/blur ratio in y (column)
        if ((nz > 1) & (size(p,1) > 1) & (size(p,2) > 1) & (size(p,3) == 1))
            areap = reshape(repmat(p.*area, 1, nz), nx, ny, nz);
        else
            areap = p .* area;      % direct multiplier of image values (p) times row ratio times column ratio
        end;
        img(ii-minx+1:ii-minx+nx, jj-miny+1:jj-miny+ny, :) = ...
            img(ii-minx+1:ii-minx+nx, jj-miny+1:jj-miny+ny, :) + (areap .* imgp) ;
    end;
end;

% Sets x and y coordinates:  make sure length(x) = number of rows in img
% and length(y) = number of cols in img
x = linspace(cx(1,1)+minx*(cx(2,1)-cx(1,1)),cx(nx,1)+maxx*(cx(2,1)-cx(1,1)),imgNrow);
y = linspace(cy(1,1)+miny*(cy(1,2)-cy(1,1)),cy(1,ny)+maxy*(cy(1,2)-cy(1,1)),imgNcol);

%% Applies Reduction (reduction in image size)
if (reduce == 1)
    nxx = length(x);
    while (x(nxx) > cx(nx,1))
        if (reduce == 1) if (sum(sum(abs(img(:,nxx,:)))) > 0) break; end; end;
        nxx = nxx - 1;
    end;
    ixx = 1;
    while (x(ixx) < cx(1,1))
        if (reduce == 1) if (sum(sum(abs(img(:,ixx,:)))) > 0) break; end; end;
        ixx = ixx + 1;
    end;
    nyy = length(y);
    while (y(nyy) > cy(1,ny))
        if (reduce == 1) if (sum(sum(abs(img(nyy,ixx:nxx,:)))) > 0) break; end; end;
        nyy = nyy - 1;
    end;
    iyy = 1;
    while (y(iyy) < cy(1,1))
        if (reduce == 1) if (sum(sum(abs(img(iyy,ixx:nxx,:)))) > 0) break; end; end;
        iyy = iyy + 1;
    end;
    x = x(ixx:nxx);
    y = y(iyy:nyy);
    img = img(ixx:nxx,iyy:nyy,:);
elseif (reduce == 2)
    img = img(1-minx:nx-minx,1-miny:ny-miny,:);
elseif (reduce == 3)
    img = img(1-minx:nx-minx,1-miny:ny-miny,:);
    img = uint8((255/double(max(img(:)))).*double(img));
end;
