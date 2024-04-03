clear; close all; clc;

addpath(genpath('./mfiles/'));

M1 = -30e-3 / 1500e-3; % front 4-f system
M2 = -100e-3 / 40e-3; % relay lens
f_MLA = 24e-3; % focal length of MLA
D_ap = 0.1524; % 6" diameter telescope aperture

pixel_pitch_ev = 15e-6; % 15 micrometer pitch event pixels

wvl = 635e-9; % wavelength for Strehl computation

pix2rad = pixel_pitch_ev .* abs(M1) ./ (abs(M2) .* f_MLA);

dataset_flag = 4;

%% load

filename_geometry = 'WFS_Geometry.mat';

% WFS geometry file
d_Geometry = load(filename_geometry);

% list the .mat files with predictions
% dd = dir('./Predictions/*.mat');
dd = dir('../Formatted_Data/*.mat');
% dd = dir('../Predictions*.mat');

names_tmp = cat(1, dd.name);
for ii = 1:size(names_tmp, 1)
    if strfind(names_tmp(ii,:), ['dataset', num2str(dataset_flag, '%02d')])
        idx_process = ii;
    end
end

% load the prediction file
dH5 = load(fullfile(dd(idx_process).folder, dd(idx_process).name), ...
    'evFrame_ref', 'flirFrame_ref', 'xc_ref_ev', 'yc_ref_ev', ...
    'xc_ref_flir', 'yc_ref_flir');
disp(dd(idx_process).name);

clear('idx_process');

% assign variables
nsub = d_Geometry.nsub;
nact = d_Geometry.nact;
xsub = d_Geometry.xsub;
ysub = d_Geometry.ysub;
xact = d_Geometry.xact;
yact = d_Geometry.yact;

% % positions_FLIR = dH5.positions;
% positions_FLIR = dH5.positions_FLIR;
% positions_baseline = dH5.positions_baseline;
% predictions_CNN = dH5.predictions_CNN;

% slope measurement indices of defined geometry
idx_slope_G = [1:nsub, 1:nsub]';

dx_act = xact(2) - xact(1);
dy_act = dx_act; % assume square subapertures

%%

evFrame_ref = dH5.evFrame_ref;
flirFrame_ref = dH5.flirFrame_ref;

xc_ref_ev = dH5.xc_ref_ev;
yc_ref_ev = dH5.yc_ref_ev;
xc_ref_flir = dH5.xc_ref_flir;
yc_ref_flir = dH5.yc_ref_flir;

figure;
imagesc(evFrame_ref);
axis image xy;
colormap(gray);
colorbar;
hold on;
% plot(dH5.xc_ref_ev, dH5.yc_ref_ev, 'rx', 'LineWidth', 2);
for ii = 1:length(xc_ref_ev)
    text(xc_ref_ev(ii), yc_ref_ev(ii), num2str(ii), 'Color', 'green', ...
        'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end
title('Event Sensor Reference Image');
% legend('Subapertures', 'location', 'northeast');
hold off;

% figure;
% imagesc(flirFrame_ref);
% axis image xy;
% colormap(gray);
% colorbar;
% hold on;
% plot(xc_ref_flir, yc_ref_flir, 'rx', 'LineWidth', 2);
% title('Frame Sensor Reference Image');
% % legend('Subapertures', 'location', 'northeast');
% hold off;

aogeom(filename_geometry);
hold on;
for ii = 1:nsub
    text(xsub(ii), ysub(ii), num2str(ii), 'Color', 'Red', 'FontSize', 12, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end
for ii = 1:length(xact)
    text(xact(ii) + 0.003, yact(ii), num2str(ii), 'Color', ...
        'White', 'FontSize', 12, 'FontWeight', 'bold');
end
xlabel('meters');
ylabel('meters');
fig1 = get(gcf);
hold off;



% exportgraphics(fig1, 'aogeom_test.pdf');

reset(groot); % reset plotting settings (they are modified by aogeom)

%% try to create custom geometry matrix G

% idx_rmv_actuators = [1:6, 7, 14, 15, 24, 56, 61, 66, 70, 71:84]; % dataset01
% idx_rmv_subapertures = [1, 5, 6, 12, 13, 20, 21, 22, 26, 27, 28, 31, 32, 33, 34, 37, 38, 39, 40, 44, 47, 48, 49, 55, 56, 60]; % dataset01

% idx_rmv_actuators = [1:6, 7, 14, 15, 24, 61, 70, 71, 78, 79:84]; % dataset02, 03, 04
% idx_rmv_subapertures = [1, 5, 6, 12, 13, 21, 22, 27, 28, 31, 33, 34, 39, 40, 44, 48, 49, 55, 56, 60]; % dataset02, 03, 04

% idx_rmv_actuators = [1:6, 7, 14, 15, 24, 25, 61, 70, 71, 78, 79:84]; % dataset07
% idx_rmv_subapertures = [1, 2, 5, 6, 12, 13, 21, 22, 27, 28, 31, 33, 34, 39, 40, 44, 48, 49, 55, 56, 58, 60]; % dataset07

% idx_rmv_actuators = [1:6, 7, 14, 15, 24, 25, 47, 56, 61, 66, 70, 71:84]; % dataset08
% idx_rmv_subapertures = [1, 2, 5, 6, 12, 13, 20, 21, 22, 26, 27, 28, 31:34, 37:40, 44, 45, 47:49, 55, 56, 60]; % dataset08

% idx_rmv_actuators = [1:6, 7, 14, 15, 24, 25, 34, 39, 47, 56, 61, 66, 70, 71:84]; % dataset10
% idx_rmv_subapertures = [1, 2, 5, 6, 12, 13, 20, 21, 22, 26, 27, 28, 31:34, 37:40, 43, 44, 45, 47, 48, 49, 55, 56, 57, 60]; % dataset10

% idx_rmv_actuators = [1:6, 7, 14, 15, 23, 24, 25, 34, 42, 50, 60, 61, 70, 71, 78, 79:84]; % dataset11
% idx_rmv_subapertures = [1, 2, 5, 6, 12, 13, 21, 22, 27, 28, 31, 33, 34, 35, 39, 40, 44, 48, 49, 50, 55, 56:60]; % dataset11

% idx_rmv_actuators = [1:6, 7, 14, 15, 24, 25, 35, 43, 51, 61, 70, 71, 78, 79:84]; % dataset13
% idx_rmv_subapertures = [1:5, 6, 12, 13, 19, 21, 22, 27, 28, 31, 33, 34, 39, 40, 44, 48, 49, 55, 56, 60]; % dataset13

% idx_rmv_actuators = [1:14, 15, 24, 61, 70, 71, 78, 79:84]; % dataset14
% idx_rmv_subapertures = [1, 5, 6, 12, 13, 14, 21, 22, 23, 27, 28, 29, 31, 33, 34, 35, 38, 39, 40, 41, 44, 48, 49, 55, 56, 60]; % dataset14

% idx_rmv_actuators = [1:14, 15, 24, 25, 34, 42, 50, 60, 61, 70, 71, 79:84]; % dataset15
% idx_rmv_subapertures = [1, 2, 5, 6, 12, 13, 14, 21:23, 27:29, 31, 33:35, 39:41, 44, 48, 49, 56:60]; % dataset15

% idx_rmv_actuators = [1:6, 7, 14, 15, 24, 25, 47, 61, 70, 71, 78, 79:84]; % dataset16
% idx_rmv_subapertures = [1, 2, 5, 6, 12, 13, 21, 22, 27, 28, 31, 33, 34, 39, 40, 44, 45, 48, 49, 55, 56, 60]; % dataset16

% idx_rmv_actuators = [1:7, 14, 15, 24, 25, 34, 42, 47, 50, 56, 60, 61, 66, 70, 71:84]; % dataset17
% idx_rmv_subapertures = [1, 2, 5, 6, 12, 13, 20, 21, 22, 26, 27, 28, 31:35, 37:40, 44, 45, 47, 48, 49, 55:60]; % dataset17

% idx_rmv_actuators = [1:6, 7, 14, 15, 24, 25, 61, 70, 71, 78, 79:84]; % dataset19
% idx_rmv_subapertures = [1, 2, 5, 6, 12, 13, 21, 22, 27, 28, 31, 33, 34, 39, 40, 44, 48, 49, 55, 56, 58, 60]; % dataset19

% idx_rmv_actuators = [1:15, 24, 25, 61, 70, 71, 78:84]; % dataset22
% idx_rmv_subapertures = [1, 2, 5, 6, 12, 13, 14, 21:23, 27:29, 31, 33:35, 39:41, 44, 48, 49, 55, 56, 60]; % dataset22

% idx_rmv_actuators = [1:7, 14, 15, 24, 25, 47, 56, 61, 66, 70:84]; % dataset24
% idx_rmv_subapertures = [1, 2, 5, 6, 12, 13, 20:22, 26:28, 31:34, 37:40, 44, 45, 47:49, 55, 56, 60]; % dataset24

% idx_rmv_actuators = [1:7, 14, 15, 24, 25, 47, 56, 61, 70, 71, 78:84]; % dataset25
% idx_rmv_subapertures = [1, 2, 5, 6, 12, 13, 21, 22, 27, 28, 31, 33, 34, 37, 39, 40, 44, 45, 48, 49, 55, 56, 58, 60]; % dataset25

% idx_rmv_actuators = [1:7, 13, 14, 15, 24, 25, 61, 70, 71, 78:84]; % dataset26
% idx_rmv_subapertures = [1, 2, 5, 6, 12, 13, 21, 22, 27, 28, 31, 33, 34, 39, 40, 41, 44, 48, 49, 55, 56, 60]; % dataset26

% idx_rmv_actuators = [1:7, 14, 15, 24, 25, 61, 70, 71, 78:84]; % dataset27
% idx_rmv_subapertures = [1, 2, 5, 6, 12, 13, 21, 22, 27, 28, 31, 33, 34, 39, 40, 44, 48, 49, 55, 56, 58, 60]; % dataset27

% idx_rmv_actuators = [1:7, 14, 15, 24, 25, 61, 70, 79:84]; % dataset28
% idx_rmv_subapertures = [1, 2, 5, 6, 13, 21, 22, 27, 28, 31, 33, 34, 39, 40, 44, 48, 49, 56, 58, 60]; % dataset28

% idx_rmv_actuators = [1:7, 14, 15, 24, 25, 61, 70, 79:84]; % dataset29
% idx_rmv_subapertures = [1, 2, 5, 6, 13, 21, 22, 27, 28, 31, 33, 34, 39, 40, 44, 48, 49, 56, 58, 60]; % dataset29

% idx_rmv_actuators = [1:15, 24, 25, 47, 56, 61, 70, 79:84]; % dataset30
% idx_rmv_subapertures = [1, 2, 5, 6, 13, 14, 21:23, 27:29, 31, 33:35, 37, 39:41, 44, 45, 48, 49, 56, 60]; % dataset30

% idx_rmv_actuators = [1:7, 14, 15, 24, 25, 61, 70, 79:84]; % dataset31
% idx_rmv_subapertures = [1, 2, 5, 6, 13, 21, 22, 27, 28, 31, 33, 34, 39, 40, 44, 48, 49, 56, 58, 60]; % dataset31

% idx_rmv_actuators = [1:7, 15, 24, 61, 70, 79:84]; % dataset32
% idx_rmv_subapertures = [1, 5, 6, 13, 21, 22, 27, 28, 31, 33, 34, 39, 40, 44, 48, 56, 58, 60]; % dataset32

% idx_rmv_actuators = [1:7, 14, 15, 24, 25, 34, 42, 47, 50, 55, 56, 60, 61, 65, 66, 70:84]; % dataset33
% idx_rmv_subapertures = [1, 2, 5, 6, 12, 13, 20:22, 25:28, 31:34, 37:40, 44, 45, 47:49, 55:60]; % dataset33

% idx_rmv_actuators = [1:15, 24, 25, 61, 70, 79:84]; % dataset34
% idx_rmv_subapertures = [1, 2, 5, 6, 13, 14, 21:23, 27:29, 31, 33:35, 39:41, 44, 48, 49, 56, 60]; % dataset34

% idx_rmv_actuators = [1:7, 14, 15, 24, 25, 61, 70, 79:84]; % dataset35
% idx_rmv_subapertures = [1, 2, 5, 6, 13, 21, 22, 27, 28, 31, 33, 34, 39, 40, 44, 48, 49, 56, 58, 60]; % dataset35


% % use this to remove top and bottom rows of subapertures and actuators
% idx_rmv_actuators = [1:6, 79:84];
% idx_rmv_subapertures = [13, 21, 22, 27, 28, 33, 34, 39, 40, 48];


xact_custom = xact; yact_custom = yact;
xsub_custom = xsub; ysub_custom = ysub;

% xact_custom(idx_rmv_actuators) = [];
% yact_custom(idx_rmv_actuators) = [];
% xsub_custom(idx_rmv_subapertures) = [];
% ysub_custom(idx_rmv_subapertures) = [];

[G, ~] = lsfptosw(xsub_custom, ysub_custom, 1, xact_custom, yact_custom);
% [G, ~] = lsfptos(xsub_custom, ysub_custom, xact_custom, yact_custom);
G = G / dx_act;


aogeom(filename_geometry);
hold on;
for ii = 1:length(xsub_custom)
    text(xsub_custom(ii), ysub_custom(ii), num2str(ii), 'Color', 'green', ...
        'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end
hold off;

aogeom(filename_geometry);
hold on;
for ii = 1:length(xact_custom)
    text(xact_custom(ii) + 0.003, yact_custom(ii), num2str(ii), 'Color', ...
        'White', 'FontSize', 12, 'FontWeight', 'bold');
end
hold off;

reset(groot); % reset plotting settings (they are modified by aogeom)

%% inspect G

figure;
imagesc(G);
axis image;
colorbar;

% figure;
% % imagesc(G(1:40, :));
% imagesc(G(1:60, :));
% axis image;
% colorbar;
% title('X-Slope');
% 
% figure;
% % imagesc(G(41:80, :));
% imagesc(G(61:120, :));
% axis image;
% colorbar;
% title('Y-Slope');
% 
% figure;
% % imagesc(G(81:end, :));
% imagesc(G(121:end, :));
% axis image;
% colorbar;
% title('Waffle');


fig1 = figure('Position', [100, 100, 1600, 400]);
subplot(1, 3, 1);
imagesc(G(1:60, :));
axis image off;
% colorbar;
title('X-Slope');

subplot(1, 3, 2);
imagesc(G(61:120, :));
axis image off;
% colorbar;
title('Y-Slope');

subplot(1, 3, 3);
imagesc(G(121:end, :));
axis image off;
% colorbar;
title('Waffle');

exportgraphics(fig1, 'TestSave.pdf');


fig1 = figure('Position', [100, 100, 1050, 400]);
subplot(1, 2, 1);
imagesc(G(1:60, :));
axis image off;
% colorbar;
title('X-Slope');

subplot(1, 2, 2);
imagesc(G(61:120, :));
axis image off;
% colorbar;
title('Y-Slope');

exportgraphics(fig1, 'TestSave.pdf');

%% compute reconstruction matrix

% compute singular values of reconstruction matrix
rsv = aorecon(G);

% plot singular values
fig1 = figure;
semilogy(rsv, '-o', 'LineWidth', 2);
xlim([0, length(rsv)+1]);
xlabel('Singular Mode Index');
ylabel('Singular Value');
set(gca, 'FontWeight', 'bold');

exportgraphics(fig1, 'Figure_FriedGeometryMatrix_WithWaffle_SingularValues.pdf');

num_singular_values_remove = 1;

% compute reconstruction matrix removing n singular modes
[Reconstructor_TiltIncluded, ~, ru, rs, rv] = aorecon(G, num_singular_values_remove, 0); % ru, rs, rv are USV from SVD
[Reconstructor_TiltRemoved, ~, ~, ~, ~] = aorecon(G, num_singular_values_remove, 1); % ru, rs, rv are same

% [u, s, v] = svd(Reconstructor_TiltIncluded);

Reconstructor_TiltIncluded = Reconstructor_TiltIncluded(:, 1:(2 * size(Reconstructor_TiltIncluded, 2) / 3));
Reconstructor_TiltRemoved = Reconstructor_TiltRemoved(:, 1:(2 * size(Reconstructor_TiltRemoved, 2) / 3));

% [u, s, v] = svd(Reconstructor_TiltIncluded);
% [u, s, v] = svd(Reconstructor_TiltRemoved);
% 
% figure;
% plot(diag(s), '-o', 'LineWidth', 2);

G_recon = G(1:(2 * size(G, 1) / 3), :); % for going back to slopes

condition_number_initial = rsv(1) / rsv(end);
condition_number_new = rsv(1) / rsv(end-num_singular_values_remove);

figure;
imagesc(Reconstructor_TiltIncluded);
axis image;
colorbar;
title('Tilt-Included Reconstructor');

figure;
imagesc(Reconstructor_TiltRemoved);
axis image;
colorbar;
title('Tilt-Removed Reconstructor');

%% saving

% clear('dH5', 'ii');
% 
% save('Reconstructor_dataset01_tmp.mat');

%% visualize geometry matrix modes

% [u, s, v] = svd(G);
% 
% figure;
% plot(diag(rs), '-*', 'LineWidth', 2);
% hold on;
% plot(diag(s), '-o', 'LineWidth', 2);


figure;
% cnt = size(rv, 2) - 20;
cnt = 0;
cnt2 = 0; % for plotting
for j = 1:5 % row
    for k = 1:4 % column
        cnt = cnt + 1;
        cnt2 = cnt2 + 1;

        % [p, x1_tmp, y1_tmp] = actmap(rv(:,cnt), xact, yact);
        % [p, x1_tmp, y1_tmp] = actmap(rwv(:,cnt), xact, yact);
        [p, x1_tmp, y1_tmp] = actmap(rv(:,cnt), xact_custom, yact_custom);
        % [p, x1_tmp, y1_tmp] = actmap(v(:,cnt), xact_custom, yact_custom);

        subplot(5, 4, cnt2);
        imagesc(x1_tmp, y1_tmp, p);
        axis image;
        colormap(jet);
        colorbar; %clim([min(rv(:)), max(rv(:))]);
        grid off;
        title(['Mode ', num2str(cnt)]);
    end
end

%% visualize reconstructor modes

[u, s, v] = svd(Reconstructor_TiltIncluded);
% [u, s, v] = svd(Reconstructor_TiltRemoved);

% figure;
% plot(diag(s), '-o', 'LineWidth', 2);


% figure;
% cnt = 0;
% % cnt = size(u, 1) - 20;
% cnt2 = 0; % for plotting
% for j = 1:5 % row
%     for k = 1:4 % column
%         cnt = cnt + 1;
%         cnt2 = cnt2 + 1;
% 
%         % [p, x1_tmp, y1_tmp] = actmap(rv(:,cnt), xact, yact);
%         % [p, x1_tmp, y1_tmp] = actmap(rwv(:,cnt), xact, yact);
%         % [p, x1_tmp, y1_tmp] = actmap(rv(:,cnt), xact_custom, yact_custom);
%         % [p, x1_tmp, y1_tmp] = actmap(v(:,cnt), xact_custom, yact_custom);
%         [p, x1_tmp, y1_tmp] = actmap(u(:,cnt), xact_custom, yact_custom);
% 
%         subplot(5, 4, cnt2);
%         imagesc(x1_tmp, y1_tmp, p);
%         axis image xy square;
%         colormap(jet);
%         colorbar; %clim([min(rv(:)), max(rv(:))]);
%         grid off;
%         title(['Mode ', num2str(cnt)]);
%     end
% end

fig1 = figure('Position', [100, 100, 1000, 700]);
cnt = 0;
% cnt = size(u, 1) - 16;
cnt2 = 0; % for plotting
for j = 1:4 % row
    for k = 1:4 % column
        cnt = cnt + 1;
        cnt2 = cnt2 + 1;

        % [p, x1_tmp, y1_tmp] = actmap(rv(:,cnt), xact, yact);
        % [p, x1_tmp, y1_tmp] = actmap(rwv(:,cnt), xact, yact);
        % [p, x1_tmp, y1_tmp] = actmap(rv(:,cnt), xact_custom, yact_custom);
        % [p, x1_tmp, y1_tmp] = actmap(v(:,cnt), xact_custom, yact_custom);
        [p, x1_tmp, y1_tmp] = actmap(u(:,cnt), xact_custom, yact_custom);

        subplot(4, 4, cnt2);
        imagesc(x1_tmp, y1_tmp, p);
        axis image xy square;
        colormap(jet);
        colorbar; %clim([min(rv(:)), max(rv(:))]);
        grid off;
        title(['Mode ', num2str(cnt)]);
    end
end
exportgraphics(fig1, 'Reconstructor_Modes.pdf');