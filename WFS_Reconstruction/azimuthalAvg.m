function [r, dataAvg] = azimuthalAvg(x, y, g)
% [r, dataAvg] = azimuthalAvg(x, y, g);
% [r, dataAvg] = azimuthalAvg(grid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the radial average of the input data field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% x [vector] = Vector of x sample locations (arb)
% y [vector] = Vector of y sample locations (arb)
% g [matrix] = 2D data to be radially averaged (arb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% r [vector] = Radial sample points (arb)
% dataAvg [vector] = Radially averaged data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% original grid info:
dx = x(2) - x(1);
dy = y(2) - y(1);
xmax = max(x);
ymax = max(y);
Nx = numel(x);
Ny = numel(y);
% make a 2D grid:
[xx, yy] = meshgrid(x, y);
% make a grid in polar coords for interpolation:
dTheta = min([atan(dx/xmax), atan(dy/ymax)]); % smallest angle
Theta = linspace(-pi, pi, 2*pi/dTheta); % range of angles
r = linspace(0, sqrt(xmax^2+ymax^2), 2*max([Nx Ny])); % range of r

% interpolate to polar grid:
xi = bsxfun(@times, r, cos(Theta).');
yi = bsxfun(@times, r, sin(Theta).');
dataPolar = interp2(xx, yy, g, xi, yi, '*linear', NaN);

dataAvg = mean(dataPolar,'omitnan');

% % find invalid values:
% nanMap = isnan(dataPolar);
% % average over valid values:
% dataAvg = nan(1, numel(r));
% for idxCol = 1 : numel(r) % loop over columns
%     % average down the rows:
%     dataAvg(idxCol) = mean(dataPolar(~nanMap(:,idxCol), idxCol));
% end
