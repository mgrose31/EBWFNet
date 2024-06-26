clear; close all; clc;

M1 = -30e-3 / 1500e-3; % front 4-f system
M2 = -100e-3 / 40e-3; % relay lens
f_MLA = 24e-3; % focal length of MLA
D_ap = 0.1524; % 6" diameter telescope aperture

pixel_pitch_ev = 15e-6; % 15 micrometer pitch event pixels

wvl = 635e-9; % wavelength for Strehl computation

pix2rad = pixel_pitch_ev .* abs(M1) ./ (abs(M2) .* f_MLA);

% dataset ID
dataset_ID = 4;

%% load reconstructor

d_reconstructor = load(['Reconstructors_DISTRIBUTE\Reconstructor_dataset', num2str(dataset_ID, '%02d'), '_DISTRIBUTE.mat']);

%% load the frame and AR(1) data

% dataset01 data
dd_data = dir(['../Data_DISTRIBUTE/Dataset', num2str(dataset_ID, '%02d'), '/*.mat']);

% pre-allocate x-axis and y-axis slopes for frame and AR1 data
sx_Frame = [];
sy_Frame = [];
sx_AR1 = [];
sy_AR1 = [];

% load the data
for ii = 1:length(dd_data)
    d_mat = load(fullfile(dd_data(ii).folder, dd_data(ii).name));
    sx_Frame(ii,:) = d_mat.xPos_Frame;
    sy_Frame(ii,:) = d_mat.yPos_Frame;
    sx_AR1(ii,:) = d_mat.xPos_AR1;
    sy_AR1(ii,:) = d_mat.yPos_AR1;
end

clear('d_mat');

% convert pixels of displacement to radians
sx_Frame = sx_Frame .* pix2rad;
sy_Frame = sy_Frame .* pix2rad;
sx_AR1 = sx_AR1 .* pix2rad;
sy_AR1 = sy_AR1 .* pix2rad;

% stack the slopes
s_Frame = [sx_Frame; sy_Frame];
s_AR1 = [sx_AR1; sy_AR1];

%% reconstruct the phase

R = d_reconstructor.Reconstructor_TiltIncluded;  % full-aptilt included
% R = d_reconstructor.Reconstructor_TiltRemoved; % full-ap tilt removed

% reconstruct the phase
phase_Frame = reconstruct_phase(R, s_Frame, wvl);
phase_AR1 = reconstruct_phase(R, s_AR1, wvl);

%% visualize the phase

% actuator positions
xact = d_reconstructor.xact_custom;
yact = d_reconstructor.yact_custom;

% convert to 2D map with padding 2
[phase_Frame_2D, x1_Frame_tmp, y1_Frame_tmp] = actmap(phase_Frame(:,1), xact, yact, 2);

% show first frame
figure;
imagesc(x1_Frame_tmp, y1_Frame_tmp, phase_Frame_2D);
axis image xy;
colormap('turbo');
colorbar;

% show a video
fig = figure;
for ii = 1:100
    [phase_Frame_2D, x1_Frame_tmp, y1_Frame_tmp] = actmap(phase_Frame(:,ii), xact, yact, 2);
    imagesc(x1_Frame_tmp, y1_Frame_tmp, phase_Frame_2D);
    axis image xy;
    colormap('turbo');
    colorbar;
    title(['Frame ', num2str(ii)]);

    drawnow;

    pause(0.1);

    clf(fig);
end

%% functions

function OPD = reconstruct_phase(R, s, wvl)
% R is reconstruction matrix
% s is slope vector (all x slopes then all y slopes)

OPD_m = R * s;  % OPD in meters
OPD = OPD_m * (2*pi) / wvl; % OPD in phase (radians)

end

function [dm, x, y] = actmap(d, xact, yact, padact)
% convert actuator positions to a 2D map with padding padact
% d is the OPD or phase

reltol = 1.0e-06;
nc = prod(size(xact));
dxs = abs(xact(2:nc)-xact(1:(nc-1)));
rdxs = reltol * max(dxs);
idxs = find(dxs > rdxs);
dxs = dxs(idxs);
dx = min(dxs);
dys = abs(yact(2:nc)-yact(1:(nc-1)));
rdys = reltol * max(dys);
idys = find(dys > rdys);
dys = dys(idys);
dy = min(dys);
eps = reltol * min([dx,dy]);
x = xact(1);
y = yact(1);
for i=2:nc
    if (length(find(abs(x-xact(i)) < eps)) <= 0)
        x = [x, xact(i)];
    end;
    if (length(find(abs(y-yact(i)) < eps)) <= 0)
        y = [y, yact(i)];
    end;
end;
x = sort(x);
y = sort(y);
nx = prod(size(x));
ny = prod(size(y));
dm = NaN(nx,ny);
for i=1:nc
    ix = find(abs(x-xact(i)) < eps);
    iy = find(abs(y-yact(i)) < eps);
    dm(ix,iy) = d(i);
end;
if (~exist('padact','var'))
    padact = 0;
end;
if (isempty(padact))
    padact = 0;
end;
for ii=1:padact
    x = [x(1)-(x(2)-x(1)), x, x(end)+(x(end)-x(end-1))];
    y = [y(1)-(y(2)-y(1)), y, y(end)+(y(end)-y(end-1))];
end;
if (padact > 0)
    dmhold = dm;
    dm = repmat(nan, nx+2*padact, ny+2*padact);
    dm(1+padact:end-padact,1+padact:end-padact) = dmhold;
end;
return
end