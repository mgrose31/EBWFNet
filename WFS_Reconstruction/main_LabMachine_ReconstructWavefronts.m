clear; close all; clc;

addpath(genpath('./mfiles/'));

M1 = -30e-3 / 1500e-3; % front 4-f system
M2 = -100e-3 / 40e-3; % relay lens
f_MLA = 24e-3; % focal length of MLA
D_ap = 0.1524; % 6" diameter telescope aperture

pixel_pitch_ev = 15e-6; % 15 micrometer pitch event pixels

wvl = 635e-9; % wavelength for Strehl computation

pix2rad = pixel_pitch_ev .* abs(M1) ./ (abs(M2) .* f_MLA);

dataset_flag = 35;

%% load

filename_geometry = 'WFS_Geometry.mat';

% WFS geometry file
d_Geometry = load(filename_geometry);

dd_Reconstructors = dir('./Reconstructors/*.mat');
recon_names_tmp = cat(1, dd_Reconstructors.name);
for ii = 1:size(recon_names_tmp, 1)
    if strfind(recon_names_tmp(ii,:), ['dataset', num2str(dataset_flag, '%02d')])
        idx_process_recon = ii;
    end
end

dH5_Reconstructor = load(fullfile(dd_Reconstructors(idx_process_recon).folder, dd_Reconstructors(idx_process_recon).name));
disp(['Loading ', dd_Reconstructors(idx_process_recon).name]);

% list the .mat files with predictions
dd_Predictions = dir('..\Trained_Networks\LargeNetwork_SlopeMSE_CartesianMeshgrid_10msIE\depth4\predictions\predictions*.mat');
pred_names_tmp = cat(1, dd_Predictions.name);
for ii = 1:size(pred_names_tmp, 1)
    if strfind(pred_names_tmp(ii,:), ['dataset', num2str(dataset_flag, '%02d')])
        idx_process_pred = ii;
    end
end

% load the prediction file
dH5_Predictions = load(fullfile(dd_Predictions(idx_process_pred).folder, dd_Predictions(idx_process_pred).name));
disp(['Loading ', dd_Predictions(idx_process_pred).name]);

% assign variables
% nsub = dH5_Reconstructor.nsub;
% nact = dH5_Reconstructor.nact;
% xsub = dH5_Reconstructor.xsub;
% ysub = dH5_Reconstructor.ysub;
% xact = dH5_Reconstructor.xact;
% yact = dH5_Reconstructor.yact;

xsub_custom = dH5_Reconstructor.xsub_custom;
ysub_custom = dH5_Reconstructor.ysub_custom;
xact_custom = dH5_Reconstructor.xact_custom;
yact_custom = dH5_Reconstructor.yact_custom;
dx_act = dH5_Reconstructor.dx_act;


positions_Frame = dH5_Predictions.positions_Frame;
predictions_AR1 = dH5_Predictions.predictions_AR1;
predictions_CNN = dH5_Predictions.predictions_CNN;

% % slope measurement indices of defined geometry
% idx_slope_G = [1:nsub, 1:nsub]';

% dx_act = xact_custom(2) - xact_custom(1);
% dy_act = dx_act; % assume square subapertures
% if abs(dx_act - 0.0164) >= 1e-10
%     error("dx_act does not equal what it should!");
% end

%%

evFrame_ref = dH5_Predictions.evFrame_ref;
flirFrame_ref = dH5_Predictions.flirFrame_ref;

xc_ref_ev = dH5_Predictions.xc_ref_ev;
yc_ref_ev = dH5_Predictions.yc_ref_ev;
xc_ref_flir = dH5_Predictions.xc_ref_flir;
yc_ref_flir = dH5_Predictions.yc_ref_flir;

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

% aogeom(filename_geometry);
% hold on;
% for ii = 1:nsub
%     text(xsub(ii), ysub(ii), num2str(ii), 'Color', 'Red', 'FontSize', 12, ...
%         'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% end
% hold off;

reset(groot); % reset plotting settings (they are modified by aogeom)

%% try to create custom geometry matrix G

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

% figure;
% imagesc(G(1:40, :));
% axis image;
% colorbar;
% title('X-Slope');
% 
% figure;
% imagesc(G(41:80, :));
% axis image;
% colorbar;
% title('Y-Slope');
% 
% figure;
% imagesc(G(81:end, :));
% axis image;
% colorbar;
% title('Waffle');

%% compute reconstruction matrix

Reconstructor_TiltIncluded = dH5_Reconstructor.Reconstructor_TiltIncluded;
Reconstructor_TiltRemoved = dH5_Reconstructor.Reconstructor_TiltRemoved;

% figure;
% imagesc(Reconstructor_TiltIncluded);
% axis image;
% colorbar;
% title('Tilt-Included Reconstructor');
% 
% figure;
% imagesc(Reconstructor_TiltRemoved);
% axis image;
% colorbar;
% title('Tilt-Removed Reconstructor');

if size(Reconstructor_TiltIncluded) ~= [length(xact_custom), length(xsub_custom)]
    error("Issue");
end

%% get time average stats

numFrames = size(positions_Frame, 3);
nsub = length(xsub_custom);

sx_FLIR_evPixels = squeeze(positions_Frame(1,:,:));
sy_FLIR_evPixels = squeeze(positions_Frame(2,:,:));
s_FLIR_evPixels = [sx_FLIR_evPixels; sy_FLIR_evPixels];

sx_Baseline_evPixels = squeeze(predictions_AR1(1,:,:));
sy_Baseline_evPixels = squeeze(predictions_AR1(2,:,:));
s_Baseline_evPixels = [sx_Baseline_evPixels; sy_Baseline_evPixels];

sx_CNN_evPixels = squeeze(predictions_CNN(1,:,:));
sy_CNN_evPixels = squeeze(predictions_CNN(2,:,:));
s_CNN_evPixels = [sx_CNN_evPixels; sy_CNN_evPixels];

s_FLIR_evPixels_timeavg = mean(s_FLIR_evPixels, 2);
s_Baseline_evPixels_timeavg = mean(s_Baseline_evPixels, 2);
s_CNN_evPixels_timeavg = mean(s_CNN_evPixels, 2);

% figure;
% plot(1:nsub, s_FLIR_evPixels_timeavg(1:nsub), '-o', 'LineWidth', 2);
% hold on;
% plot(1:nsub, s_Baseline_evPixels_timeavg(1:nsub), '-o', 'LineWidth', 2);
% plot(1:nsub, s_CNN_evPixels_timeavg(1:nsub), '-o', 'LineWidth', 2);
% xlabel('Subaperture');
% ylabel('x-slope time average (event pixels)');
% legend('FLIR', 'Baseline', 'CNN', 'location', 'best');
% hold off;
% 
% figure;
% plot(1:nsub, s_FLIR_evPixels_timeavg(nsub+1:end), '-o', 'LineWidth', 2);
% hold on;
% plot(1:nsub, s_Baseline_evPixels_timeavg(nsub+1:end), '-o', 'LineWidth', 2);
% plot(1:nsub, s_CNN_evPixels_timeavg(nsub+1:end), '-o', 'LineWidth', 2);
% xlabel('Subaperture');
% ylabel('y-slope time average (event pixels)');
% legend('FLIR', 'Baseline', 'CNN', 'location', 'best');
% hold off;

flag_rmv_timeavg = true;
% flag_rmv_timeavg = false;
if flag_rmv_timeavg
    disp("Removing time-average slopes.");

    s_FLIR_evPixels = s_FLIR_evPixels - repmat(s_FLIR_evPixels_timeavg, [1, numFrames]);
    s_Baseline_evPixels = s_Baseline_evPixels - repmat(s_Baseline_evPixels_timeavg, [1, numFrames]);
    s_CNN_evPixels = s_CNN_evPixels - repmat(s_CNN_evPixels_timeavg, [1, numFrames]);
end

%% get per-frame stats (full-aperture tilt)

sx_FLIR_fullApTilt = mean(s_FLIR_evPixels(1:nsub, :), 1);
sy_FLIR_fullApTilt = mean(s_FLIR_evPixels(nsub+1:end, :), 1);

sx_Baseline_fullApTilt = mean(s_Baseline_evPixels(1:nsub, :), 1);
sy_Baseline_fullApTilt = mean(s_Baseline_evPixels(nsub+1:end, :), 1);

sx_CNN_fullApTilt = mean(s_CNN_evPixels(1:nsub, :), 1);
sy_CNN_fullApTilt = mean(s_CNN_evPixels(nsub+1:end, :), 1);

% figure('Position', [100, 100, 1000, 400]);
% plot(1:numFrames, sx_FLIR_fullApTilt, '-', 'LineWidth', 2);
% hold on;
% plot(1:numFrames, sx_Baseline_fullApTilt, '-', 'LineWidth', 2);
% plot(1:numFrames, sx_CNN_fullApTilt, '-', 'LineWidth', 2);
% xlabel('Frame');
% ylabel('x-axis full-aperture tilt (event pixels)');
% legend('FLIR', 'Baseline', 'CNN', 'location', 'best');
% hold off;
% 
% figure('Position', [100, 100, 1000, 400]);
% plot(1:numFrames, sy_FLIR_fullApTilt, '-', 'LineWidth', 2);
% hold on;
% plot(1:numFrames, sy_Baseline_fullApTilt, '-', 'LineWidth', 2);
% plot(1:numFrames, sy_CNN_fullApTilt, '-', 'LineWidth', 2);
% xlabel('Frame');
% ylabel('y-axis full-aperture tilt (event pixels)');
% legend('FLIR', 'Baseline', 'CNN', 'location', 'best');
% hold off;

flag_rmv_fullApTilt = true;
% flag_rmv_fullApTilt = false;
if flag_rmv_fullApTilt
    disp("Manually removing full-aperture tilt.");

    s_FLIR_evPixels(1:nsub, :) = s_FLIR_evPixels(1:nsub, :) - repmat(sx_FLIR_fullApTilt, [nsub, 1]);
    s_FLIR_evPixels(nsub+1:end, :) = s_FLIR_evPixels(nsub+1:end, :) - repmat(sy_FLIR_fullApTilt, [nsub, 1]);

    s_Baseline_evPixels(1:nsub, :) = s_Baseline_evPixels(1:nsub, :) - repmat(sx_Baseline_fullApTilt, [nsub, 1]);
    s_Baseline_evPixels(nsub+1:end, :) = s_Baseline_evPixels(nsub+1:end, :) - repmat(sy_Baseline_fullApTilt, [nsub, 1]);

    s_CNN_evPixels(1:nsub, :) = s_CNN_evPixels(1:nsub, :) - repmat(sx_CNN_fullApTilt, [nsub, 1]);
    s_CNN_evPixels(nsub+1:end, :) = s_CNN_evPixels(nsub+1:end, :) - repmat(sy_CNN_fullApTilt, [nsub, 1]);

    % flag_tilt_remove = false;
else
    flag_tilt_remove = true;
    % disp("")
end

%% convert slopes to radians

s_FLIR_radians = s_FLIR_evPixels .* pix2rad;
s_Baseline_radians = s_Baseline_evPixels .* pix2rad;
s_CNN_radians = s_CNN_evPixels .* pix2rad;

%%

% flag_tilt_remove = true;
flag_tilt_remove = false;
if flag_tilt_remove
    disp("Using tilt-removed reconstructor.");
    Recon = Reconstructor_TiltRemoved;
else
    disp("Using tilt-included reconstructor.");
    Recon = Reconstructor_TiltIncluded;
end

for idx_frame = 1:numFrames

    s_FLIR_tmp = s_FLIR_radians(:, idx_frame);
    s_Baseline_tmp = s_Baseline_radians(:, idx_frame);
    s_CNN_tmp = s_CNN_radians(:, idx_frame);

    OPD_FLIR_tmp = reconstruct_phase(Recon, s_FLIR_tmp, wvl);
    OPD_baseline_tmp = reconstruct_phase(Recon, s_Baseline_tmp, wvl);
    OPD_CNN_tmp = reconstruct_phase(Recon, s_CNN_tmp, wvl);

    [Strehl_baseline_tmp, MSE_baseline_tmp] = Strehl_approximation(OPD_FLIR_tmp, OPD_baseline_tmp);
    [Strehl_CNN_tmp, MSE_CNN_tmp] = Strehl_approximation(OPD_FLIR_tmp, OPD_CNN_tmp);

    s_FLIR(:, idx_frame) = s_FLIR_tmp;
    s_baseline(:, idx_frame) = s_Baseline_tmp;
    s_CNN(:, idx_frame) = s_CNN_tmp;

    OPD_FLIR(:, idx_frame) = OPD_FLIR_tmp;
    OPD_baseline(:, idx_frame) = OPD_baseline_tmp;
    OPD_CNN(:, idx_frame) = OPD_CNN_tmp;

    MSE_baseline(1, idx_frame) = MSE_baseline_tmp;
    MSE_CNN(1, idx_frame) = MSE_CNN_tmp;

    Strehl_baseline(1, idx_frame) = Strehl_baseline_tmp;
    Strehl_CNN(1, idx_frame) = Strehl_CNN_tmp;

end

clear('OPD_FLIR_tmp', 'OPD_baseline_tmp', 'OPD_CNN_tmp', ...
    'sx_Baseline_evPixels', 'sy_baseline_evPixels', 'sx_FLIR_evPixels', ...
    'sy_FLIR_evPixels', 'sx_CNN_evPixels', 'sy_CNN_evPixels', ...
    'sx_CNN_radians', 'sy_CNN_radians', 'sx_FLIR_radians', ...
    'sy_FLIR_radians', 'sx_baseline_radians', 'sy_baseline_radians', ...
    'MSE_baseline_tmp', 'MSE_CNN_tmp', 's_FLIR_tmp', 's_Baseline_tmp', ...
    's_CNN_tmp', 'Strehl_baseline_tmp', 'Strehl_CNN_tmp');

%% Plot Strehl and wavefront MSE

% figure('Position', [100, 100, 1000, 400]);
% plot(1:numFrames, Strehl_baseline, '-o');
% hold on;
% plot(1:numFrames, Strehl_CNN, '-o');
% xlim([0, numFrames]);
% ylim([0, 1]);
% xlabel('Frame');
% ylabel('Approximated Strehl Ratio');
% legend('Baseline', 'CNN', 'location', 'southeast');
% hold off;
% 
% figure('Position', [100, 100, 1000, 400]);
% plot(1:numFrames, MSE_baseline, '-o');
% hold on;
% plot(1:numFrames, MSE_CNN, '-o');
% xlim([0, numFrames]);
% xlabel('Frame');
% ylabel('Wavefront MSE [rad^{2}]')
% legend('Baseline', 'CNN', 'location', 'northeast');
% hold off;

figure;
histogram(Strehl_baseline, 0:0.01:1);
hold on;
histogram(Strehl_CNN, 0:0.01:1);
xlim([0, 1]);
xlabel('Strehl Ratio');
ylabel('Counts');
legend('Baseline', 'CNN', 'location', 'northwest');
hold off;

%% save file for later processing

save(['dataset', num2str(dataset_flag, '%02d'), '_Reconstructed_Wavefronts.mat']);

% %% go back to slopes
% 
% G_recon = dH5_Reconstructor.G_recon;
% 
% s_FLIR_recon = G_recon * (OPD_FLIR * wvl / (2 * pi));
% s_baseline_recon = G_recon * (OPD_baseline * wvl / (2 * pi));
% s_CNN_recon = G_recon * (OPD_CNN * wvl / (2 * pi));
% 
% num_slopes = size(s_FLIR, 1) / 2;
% 
% SE_slope_recon_FLIR = (s_FLIR - s_FLIR_recon).^2;
% RSE_slope_recon_FLIR = sqrt(SE_slope_recon_FLIR(1:num_slopes,:) + SE_slope_recon_FLIR(num_slopes+1:end,:));
% 
% SE_slope_recon_baseline = (s_FLIR - s_baseline_recon).^2;
% RSE_slope_recon_baseline = sqrt(SE_slope_recon_baseline(1:num_slopes,:) + SE_slope_recon_baseline(num_slopes+1:end,:));
% 
% SE_slope_recon_CNN = (s_FLIR - s_CNN_recon).^2;
% RSE_slope_recon_CNN = sqrt(SE_slope_recon_CNN(1:num_slopes,:) + SE_slope_recon_CNN(num_slopes+1:end,:));
% 
% % figure;
% % plot(1:size(RSE_slope_recon_FLIR, 1), mean(RSE_slope_recon_FLIR, 2) .* 1e6, '-o', 'LineWidth', 2);
% % hold on;
% % plot(1:size(RSE_slope_recon_baseline, 1), mean(RSE_slope_recon_baseline, 2) .* 1e6, '-o', 'LineWidth', 2);
% % plot(1:size(RSE_slope_recon_CNN, 1), mean(RSE_slope_recon_CNN, 2) .* 1e6, '-o', 'LineWidth', 2);
% % legend('FLIR', 'Baseline', 'CNN');
% % hold off;
% % 
% % figure('Position', [100, 100, 1000, 400]);
% % plot(s_FLIR(8,:), 'LineWidth', 2);
% % hold on;
% % plot(s_FLIR_recon(8,:), 'LineWidth', 2);
% % plot(s_baseline_recon(8,:), 'LineWidth', 2);
% % plot(s_CNN_recon(8,:), 'LineWidth', 2);
% % legend('Original FLIR', 'Reconstructed FLIR', 'Baseline', 'CNN', 'location', 'best');
% % hold off;
% 
% %% compare OPD
% 
% % actuator_plt = 20;
% % 
% % figure('Position', [100, 100, 1000, 400]);
% % plot(1:numFrames, OPD_FLIR(actuator_plt,:), 'LineWidth', 2);
% % hold on;
% % plot(1:numFrames, OPD_baseline(actuator_plt,:), 'LineWidth', 2);
% % plot(1:numFrames, OPD_CNN(actuator_plt,:), 'LineWidth', 2);
% % xlim([0, numFrames+1]);
% % xlabel('Frame');
% % ylabel('Phase \phi [rad]');
% % legend('FLIR', 'Baseline', 'CNN (Ours)');
% % hold off;
% 
% SE_OPD_FLIR_baseline = (OPD_FLIR - OPD_baseline).^2;
% SE_OPD_FLIR_CNN = (OPD_FLIR - OPD_CNN).^2;
% 
% RMSE_OPD_FLIR_baseline = sqrt(mean(SE_OPD_FLIR_baseline, 'all'));
% RMSE_OPD_FLIR_CNN = sqrt(mean(SE_OPD_FLIR_CNN, 'all'));
% 
% % figure;
% % histogram(SE_OPD_FLIR_baseline, 0:0.1:10);
% % hold on;
% % histogram(SE_OPD_FLIR_CNN, 0:0.1:10);
% % xlabel('OPD Squared Error');
% % ylabel('Counts');
% % legend('Baseline', 'CNN', 'location', 'northeast');
% % hold off;
% 
% var_wavefront_FLIR = var(OPD_FLIR, [], 'all');
% var_wavefront_baseline = var(OPD_baseline, [], 'all');
% var_wavefront_CNN = var(OPD_CNN, [], 'all');
% 
% fprintf(1, ['FLIR wavefront variance ', num2str(var_wavefront_FLIR), '\n']);
% fprintf(1, ['Baseline wavefront variance ', num2str(var_wavefront_baseline), '\n']);
% fprintf(1, ['CNN wavefront variance ', num2str(var_wavefront_CNN), '\n']);
% 
% %% compute per-frame structure function
% 
% [opdrec2d_FLIR_tmp, x1_FLIR_tmp, y1_FLIR_tmp] = actmap(OPD_FLIR(:,1), xact_custom, yact_custom, [], [], 20);
% 
% mask = ~isnan(opdrec2d_FLIR_tmp);
% 
% delta = dx_act;
% 
% for ii=1:numFrames
%     [opdrec2d_FLIR_tmp, ~, ~] = actmap(OPD_FLIR(:,ii), xact_custom, yact_custom, [], [], 20);
%     [opdrec2d_Baseline_tmp, ~, ~] = actmap(OPD_baseline(:,ii), xact_custom, yact_custom, [], [], 20);
%     [opdrec2d_CNN_tmp, ~, ~] = actmap(OPD_CNN(:,ii), xact_custom, yact_custom, [], [], 20);
% 
%     opdrec2d_FLIR(:,:,ii) = opdrec2d_FLIR_tmp;
%     opdrec2d_Baseline(:,:,ii) = opdrec2d_Baseline_tmp;
%     opdrec2d_CNN(:,:,ii) = opdrec2d_CNN_tmp;
% 
%     opdrec2d_FLIR_tmp(~mask) = 0;
%     opdrec2d_Baseline_tmp(~mask) = 0;
%     opdrec2d_CNN_tmp(~mask) = 0;
% 
%     [D_FLIR_tmp, ~] = str_fcn2_ft_new(opdrec2d_FLIR_tmp, mask, delta);
%     [D_Baseline_tmp, ~] = str_fcn2_ft_new(opdrec2d_Baseline_tmp, mask, delta);
%     [D_CNN_tmp, ~] = str_fcn2_ft_new(opdrec2d_CNN_tmp, mask, delta);
% 
%     D_FLIR(:,:,ii) = D_FLIR_tmp;
%     D_Baseline(:,:,ii) = D_Baseline_tmp;
%     D_CNN(:,:,ii) = D_CNN_tmp;
% end
% 
% clear('D_FLIR_tmp', 'D_Baseline_tmp', 'D_CNN_tmp', 'opdrec2d_FLIR_tmp', ...
%     'opdrec2d_Baseline_tmp', 'opdrec2d_CNN_tmp');
% 
% %% analysis of reconstructed 2D OPD function
% 
% SE_baseline_opdrec2d = (opdrec2d_FLIR - opdrec2d_Baseline).^2;
% SE_CNN_opdrec2d = (opdrec2d_FLIR - opdrec2d_CNN).^2;
% 
% % figure;
% % histogram(SE_baseline_opdrec2d, 0:0.1:30);
% % hold on;
% % histogram(SE_CNN_opdrec2d, 0:0.1:30);
% % xlabel('OPD Squared Error');
% % ylabel('Counts');
% % legend('Baseline', 'CNN', 'location', 'northeast');
% % hold off;
% 
% SE_baseline_opdrec2d_mean = mean(SE_baseline_opdrec2d, 3);
% SE_CNN_opdrec2d_mean = mean(SE_CNN_opdrec2d, 3);
% 
% cmin_tmp = min([SE_baseline_opdrec2d_mean(:); SE_CNN_opdrec2d_mean(:)]);
% cmax_tmp = max([SE_baseline_opdrec2d_mean(:); SE_CNN_opdrec2d_mean(:)]);
% 
% 
% figure('Position', [100, 100, 1000, 400]);
% subplot(1, 2, 1);
% imagesc(x1_FLIR_tmp, y1_FLIR_tmp, SE_baseline_opdrec2d_mean);
% axis image;
% colormap(jet);
% colorbar; clim([cmin_tmp, cmax_tmp]);
% title('Baseline OPD MSE');
% 
% subplot(1, 2, 2);
% imagesc(x1_FLIR_tmp, y1_FLIR_tmp, SE_CNN_opdrec2d_mean);
% axis image;
% colormap(jet);
% colorbar; clim([cmin_tmp, cmax_tmp]);
% title('CNN OPD MSE');
% 
% 
% % figure('Position', [100, 100, 1200, 400]);
% % subplot(1, 3, 1);
% % imagesc(x1_FLIR_tmp, y1_FLIR_tmp, opdrec2d_FLIR(:,:,1442));
% % axis image;
% % colormap(jet);
% % colorbar;
% % 
% % subplot(1, 3, 2);
% % imagesc(x1_FLIR_tmp, y1_FLIR_tmp, opdrec2d_CNN(:,:,1442));
% % axis image;
% % colormap(jet);
% % colorbar;
% % 
% % subplot(1, 3, 3);
% % imagesc(x1_FLIR_tmp, y1_FLIR_tmp, SE_CNN_opdrec2d(:,:,1442));
% % axis image;
% % colormap(jet);
% % colorbar;
% 
% %% analysis of structure function
% 
% SE_D_Baseline = (D_FLIR - D_Baseline).^2;
% SE_D_CNN = (D_FLIR - D_CNN).^2;
% 
% SE_D_Baseline_mean = mean(SE_D_Baseline, 3);
% SE_D_CNN_mean = mean(SE_D_CNN, 3);
% 
% figure;
% histogram(SE_D_Baseline);
% hold on;
% histogram(SE_D_CNN);
% xlabel('Structure Function Squared Error');
% ylabel('Counts');
% legend('Baseline', 'CNN');
% 
% 
% cmin_tmp = min([SE_D_Baseline_mean(:); SE_D_CNN_mean(:)]);
% cmax_tmp = max([SE_D_Baseline_mean(:); SE_D_CNN_mean(:)]);
% 
% 
% figure('Position', [100, 100, 1000, 400]);
% subplot(1, 2, 1);
% imagesc(x1_FLIR_tmp, y1_FLIR_tmp, SE_D_Baseline_mean);
% axis image;
% colormap(jet);
% colorbar; clim([cmin_tmp, cmax_tmp]);
% title('Baseline Structure Function MSE');
% 
% subplot(1, 2, 2);
% imagesc(x1_FLIR_tmp, y1_FLIR_tmp, SE_D_CNN_mean);
% axis image;
% colormap(jet);
% colorbar; clim([cmin_tmp, cmax_tmp]);
% title('CNN Structure Function MSE');
% 
% %%
% 
% D_FLIR_avg = mean(D_FLIR, 3);
% D_Baseline_avg = mean(D_Baseline, 3);
% D_CNN_avg = mean(D_CNN, 3);
% 
% D_FLIR_var = var(D_FLIR, [], 3) ./ numFrames;
% D_Baseline_var = var(D_Baseline, [], 3) ./ numFrames;
% D_CNN_var = var(D_CNN, [], 3) ./ numFrames;
% 
% cmax_tmp = max([D_FLIR_avg(:); D_Baseline_avg(:); D_CNN_avg(:)]);
% 
% figure('Position', [100, 100, 1500, 400]);
% subplot(1, 3, 1);
% imagesc(x1_FLIR_tmp, y1_FLIR_tmp, D_FLIR_avg);
% axis image;
% colormap(jet);
% colorbar; clim([0, cmax_tmp]);
% title('FLIR Structure Function Mean');
% 
% subplot(1, 3, 2);
% imagesc(x1_FLIR_tmp, y1_FLIR_tmp, D_Baseline_avg);
% axis image;
% colormap(jet);
% colorbar; clim([0, cmax_tmp]);
% title('Baseline Structure Function Mean');
% 
% subplot(1, 3, 3);
% imagesc(x1_FLIR_tmp, y1_FLIR_tmp, D_CNN_avg);
% axis image;
% colormap(jet);
% colorbar; clim([0, cmax_tmp]);
% title('CNN Structure Function Mean');
% 
% 
% cmax_tmp = max([D_FLIR_var(:); D_Baseline_var(:); D_CNN_var(:)]);
% 
% figure('Position', [100, 100, 1500, 400]);
% subplot(1, 3, 1);
% imagesc(x1_FLIR_tmp, y1_FLIR_tmp, D_FLIR_var);
% axis image;
% colormap(jet);
% colorbar; clim([0, cmax_tmp]);
% title('FLIR Structure Function Variance');
% 
% subplot(1, 3, 2);
% imagesc(x1_FLIR_tmp, y1_FLIR_tmp, D_Baseline_var);
% axis image;
% colormap(jet);
% colorbar; clim([0, cmax_tmp]);
% title('Baseline Structure Function Variance');
% 
% subplot(1, 3, 3);
% imagesc(x1_FLIR_tmp, y1_FLIR_tmp, D_CNN_var);
% axis image;
% colormap(jet);
% colorbar; clim([0, cmax_tmp]);
% title('CNN Structure Function Variance');
% 
% % SE_baseline = (D_truth(:) - D_baseline(:)).^2;
% % SE_CNN = (D_truth(:) - D_CNN(:)).^2;
% % 
% % idx_use_metric = SE_baseline ~= 0 & SE_CNN ~= 0;
% % % idx_use_metric = idx_use_metric & (SE_baseline < 100 & SE_CNN < 100);
% % 
% % SE_baseline = SE_baseline(idx_use_metric);
% % SE_CNN = SE_CNN(idx_use_metric);
% 
% 
% SE_D_Baseline_avg = (D_FLIR_avg - D_Baseline_avg).^2;
% SE_D_CNN_avg = (D_FLIR_avg - D_CNN_avg).^2;
% 
% cmax_tmp = max([SE_D_Baseline_avg(:); SE_D_CNN_avg(:)]);
% 
% figure('Position', [100, 100, 1200, 400]);
% subplot(1, 2, 1);
% imagesc(x1_FLIR_tmp, y1_FLIR_tmp, SE_D_Baseline_avg);
% axis image;
% colormap(jet);
% colorbar;
% clim([0, cmax_tmp]);
% % clim([0, 0.1]);
% title('Baseline Average Structure Function MSE');
% 
% subplot(1, 2, 2);
% imagesc(x1_FLIR_tmp, y1_FLIR_tmp, SE_D_CNN_avg);
% axis image;
% colormap(jet);
% colorbar;
% clim([0, cmax_tmp]);
% % clim([0, 0.1]);
% title('CNN Average Structure Function MSE');
% % 
% % figure;
% % imagesc(x1_truth_tmp, y1_truth_tmp, D_baseline_avg - D_CNN_avg);
% % axis image;
% % colormap(jet);
% % colorbar;
% % title('Baseline - CNN');
% 
% %% Plot D(r)
% 
% [r_FLIR_tmp, dataAvg_FLIR_tmp] = azimuthalAvg(x1_FLIR_tmp, y1_FLIR_tmp, D_FLIR_avg');
% [r_Baseline_tmp, dataAvg_Baseline_tmp] = azimuthalAvg(x1_FLIR_tmp, y1_FLIR_tmp, D_Baseline_avg');
% [r_CNN_tmp, dataAvg_CNN_tmp] = azimuthalAvg(x1_FLIR_tmp, y1_FLIR_tmp, D_CNN_avg');
% 
% [r2_FLIR_tmp, dataAvg2_FLIR_tmp] = azimuthalAvg(x1_FLIR_tmp, y1_FLIR_tmp, D_FLIR_var');
% [r2_Baseline_tmp, dataAvg2_Baseline_tmp] = azimuthalAvg(x1_FLIR_tmp, y1_FLIR_tmp, D_Baseline_var');
% [r2_CNN_tmp, dataAvg2_CNN_tmp] = azimuthalAvg(x1_FLIR_tmp, y1_FLIR_tmp, D_CNN_var');
% 
% % r0_theory_tmp = 0.029;
% r0_theory_tmp = 0.028;
% 
% if flag_tilt_remove
%     alpha = 1; % no scintillation
%     % alpha = 0.5; % scintyillation
%     D_theory = 6.88 .* (r_FLIR_tmp ./ r0_theory_tmp).^(5/3) .* (1 - (alpha .* (r_FLIR_tmp ./ D_ap).^(1/3)));
% else
%     D_theory = 6.88 .* (r_FLIR_tmp ./ r0_theory_tmp).^(5/3);
% end
% 
% alpha = 1; % no scintillation
% % alpha = 0.5; % scintyillation
% D_theory = 6.88 .* (r_FLIR_tmp ./ r0_theory_tmp).^(5/3) .* (1 - (alpha .* (r_FLIR_tmp ./ D_ap).^(1/3)));
% 
% figure;
% errorbar(r_FLIR_tmp, dataAvg_FLIR_tmp, dataAvg2_FLIR_tmp, 'LineWidth', 2);
% hold on;
% errorbar(r_Baseline_tmp, dataAvg_Baseline_tmp, dataAvg2_Baseline_tmp, 'LineWidth', 2);
% errorbar(r_CNN_tmp, dataAvg_CNN_tmp, dataAvg2_CNN_tmp, 'LineWidth', 2);
% plot(r_FLIR_tmp, D_theory, 'LineWidth', 2);
% set(gca, 'YScale', 'log', 'XScale', 'log');
% xlim([0.01, 0.16]);
% xlabel('r [m]');
% ylabel('Structure Function D(r)');
% legend('FLIR', 'Baseline', 'CNN', 'r_{0}', 'location', 'best');
% hold off;
% 
% %% D(r) for individual wavefronts
% 
% for ii = 1:numFrames
%     [r_FLIR_tmp2, dataAvg_FLIR_tmp2(:,ii)] = azimuthalAvg(x1_FLIR_tmp, y1_FLIR_tmp, D_FLIR(:,:,ii)');
%     [r_Baseline_tmp2, dataAvg_Baseline_tmp2(:,ii)] = azimuthalAvg(x1_FLIR_tmp, y1_FLIR_tmp, D_Baseline(:,:,ii)');
%     [r_CNN_tmp2, dataAvg_CNN_tmp2(:,ii)] = azimuthalAvg(x1_FLIR_tmp, y1_FLIR_tmp, D_CNN(:,:,ii)');
% end
% 
% dataAvg_FLIR_tmp3 = mean(dataAvg_FLIR_tmp2(:,1:600), 2);
% dataAvg_Baseline_tmp3 = mean(dataAvg_Baseline_tmp2(:,1:600), 2);
% dataAvg_CNN_tmp3 = mean(dataAvg_CNN_tmp2(:,1:600), 2);
% 
% %%
% 
% r0_theory_tmp = 0.02;
% 
% % if flag_tilt_remove
% %     alpha = 1; % no scintillation
% %     % alpha = 0.5; % scintillation
% %     D_theory = 6.88 .* (r_FLIR_tmp ./ r0_theory_tmp).^(5/3) .* (1 - (alpha .* (r_FLIR_tmp ./ D_ap).^(1/3)));
% % else
% %     D_theory = 6.88 .* (r_FLIR_tmp ./ r0_theory_tmp).^(5/3);
% % end
% 
% alpha = 1; % no scintillation
% % alpha = 0.5; % scintillation
% D_theory = 6.88 .* (r_FLIR_tmp ./ r0_theory_tmp).^(5/3) .* (1 - (alpha .* (r_FLIR_tmp ./ D_ap).^(1/3)));
% 
% figure;
% plot(r_FLIR_tmp, dataAvg_FLIR_tmp3, 'LineWidth', 2);
% hold on;
% plot(r_Baseline_tmp, dataAvg_Baseline_tmp3, 'LineWidth', 2);
% plot(r_CNN_tmp, dataAvg_CNN_tmp3, 'LineWidth', 2);
% plot(r_FLIR_tmp, D_theory, 'LineWidth', 2);
% % set(gca, 'YScale', 'log', 'XScale', 'log');
% xlim([0.01, 0.16]);
% xlabel('r [m]');
% ylabel('Structure Function D(r)');
% legend('FLIR', 'Baseline', 'CNN', 'r_{0}', 'location', 'best');
% hold off;
% 
% %% plot reconstructed wavefronts
% 
% flag_mkVideo = true;
% if flag_mkVideo
%     v = VideoWriter('ReconstructedWavefronts_Video', 'MPEG-4');
%     v.Quality = 100;
%     v.FrameRate = 10;
% 
%     open(v);
% end
% 
% fig1 = figure('Position', [100, 100, 2000, 600]);
% hold on;
% 
% % cmin = min([OPD_FLIR(:); OPD_baseline(:); OPD_CNN(:)]);
% % cmax = max([OPD_FLIR(:); OPD_baseline(:); OPD_CNN(:)]);
% 
% cmin = -15;
% cmax = 15;
% 
% % for ii = 1:numFrames
% for ii = 1:500
% % for ii = 88
% 
%     [opdrec2d_FLIR, x1_FLIR, y1_FLIR] = actmap(OPD_FLIR(:,ii), xact_custom, yact_custom, [], [], 3);
%     [opdrec2d_Baseline, x1_baseline, y1_baseline] = actmap(OPD_baseline(:,ii), xact_custom, yact_custom, [], [], 3);
%     [opdrec2d_CNN, x1_NN, y1_NN] = actmap(OPD_CNN(:,ii), xact_custom, yact_custom, [], [], 3);
% 
%     % cmin = min([OPD_FLIR(:,ii); OPD_baseline(:,ii); OPD_CNN(:,ii)], [], 'all');
%     % cmax = max([OPD_FLIR(:,ii); OPD_baseline(:,ii); OPD_CNN(:,ii)], [], 'all');
% 
%     subplot(1, 3, 1);
%     imagesc(y1_FLIR, x1_FLIR, opdrec2d_FLIR'); % transpose for correct orientation
%     axis image xy square;
%     colormap(jet);
%     colorbar;
%     clim([cmin, cmax]);
%     title('FLIR');
%     xlabel('X-Axis');
%     ylabel('Y-Axis');
%     set(gca, 'FontWeight', 'bold');
% 
%     subplot(1, 3, 2);
%     imagesc(y1_baseline, x1_baseline, opdrec2d_Baseline'); % transpose for correct orientation
%     axis image xy square;
%     colormap(jet);
%     colorbar;
%     clim([cmin, cmax]);
%     title('Baseline');
%     xlabel('X-Axis');
%     ylabel('Y-Axis');
%     set(gca, 'FontWeight', 'bold');
% 
%     subplot(1, 3, 3);
%     imagesc(y1_NN, x1_NN, opdrec2d_CNN'); % transpose for correct orientation
%     axis image xy square;
%     colormap(jet);
%     colorbar;
%     clim([cmin, cmax]);
%     title('CNN (Ours)');
%     xlabel('X-Axis');
%     ylabel('Y-Axis');
%     set(gca, 'FontWeight', 'bold');
% 
%     sgtitle(['Frame ', num2str(ii)]);
% 
%     drawnow;
% 
%     if flag_mkVideo
%         writeVideo(v, getframe(fig1));
%     end
% 
%     clf(fig1);
% end
% 
% if flag_mkVideo
%     close(v);
% end

% fig1 = figure;
% imagesc(y1_FLIR, x1_FLIR, opdrec2d_FLIR'); % transpose for correct orientation
% axis image xy square;
% colormap(jet);
% colorbar;
% clim([cmin, cmax]);
% title('FLIR');
% xlabel('X-Axis');
% ylabel('Y-Axis');
% set(gca, 'FontWeight', 'bold');
% 
% fig2 = figure;
% imagesc(y1_baseline, x1_baseline, opdrec2d_Baseline'); % transpose for correct orientation
% axis image xy square;
% colormap(jet);
% colorbar;
% clim([cmin, cmax]);
% title('Baseline');
% xlabel('X-Axis');
% ylabel('Y-Axis');
% set(gca, 'FontWeight', 'bold');
% 
% fig3 = figure;
% imagesc(y1_NN, x1_NN, opdrec2d_CNN'); % transpose for correct orientation
% axis image xy square;
% colormap(jet);
% colorbar;
% clim([cmin, cmax]);
% title('CNN (Ours)');
% xlabel('X-Axis');
% ylabel('Y-Axis');
% set(gca, 'FontWeight', 'bold');
% 
% exportgraphics(fig1, 'Figure_ReconstructedWavefront_FLIR_dataset04_Frame88.pdf');
% exportgraphics(fig2, 'Figure_ReconstructedWavefront_Baseline_dataset04_Frame88.pdf');
% exportgraphics(fig3, 'Figure_ReconstructedWavefront_CNN_dataset04_Frame88.pdf');




% %%
% 
% fig1 = figure;
% hold on;
% 
% cmin = min(OPD_FLIR(:));
% cmax = max(OPD_FLIR(:));
% 
% v = VideoWriter('ReconstructedWavefront_example', 'MPEG-4');
% v.Quality = 100;
% v.FrameRate = 10;
% 
% open(v);
% 
% % for ii = 1:numFrames
% for ii = 1:600
%     [opdrec2d_FLIR, x1_FLIR, y1_FLIR] = actmap(OPD_FLIR(:,ii), xact, yact);
%     % [opdrec2d_baseline, x1_baseline, y1_baseline] = actmap(OPD_baseline(:,ii), xact, yact);
%     % [opdrec2d_CNN, x1_NN, y1_NN] = actmap(OPD_CNN(:,ii), xact, yact);
% 
%     % cmin = min([OPD_truth(:,ii); OPD_baseline(:,ii); OPD_CNN(:,ii)], [], 'all');
%     % cmax = max([OPD_truth(:,ii); OPD_baseline(:,ii); OPD_CNN(:,ii)], [], 'all');
% 
%     % subplot(1, 3, 1);
%     imagesc(x1_FLIR, y1_FLIR, opdrec2d_FLIR);
%     axis image;
%     colormap(jet);
%     colorbar;
%     clim([cmin, cmax]);
%     title(['Frame ', num2str(ii)]);
% 
%     % subplot(1, 3, 2);
%     % imagesc(x1_baseline, y1_baseline, opdrec2d_baseline);
%     % axis image;
%     % colormap(jet);
%     % colorbar;
%     % clim([cmin, cmax]);
%     % title('Baseline');
%     % 
%     % subplot(1, 3, 3);
%     % imagesc(x1_NN, y1_NN, opdrec2d_CNN);
%     % axis image;
%     % colormap(jet);
%     % colorbar;
%     % clim([cmin, cmax]);
%     % title('CNN (Ours)');
%     % 
%     % sgtitle(['Frame ', num2str(ii)]);
% 
%     drawnow;
%     writeVideo(v, getframe(fig1));
%     clf(fig1);
% end
% 
% close(v);


%% functions

function OPD = reconstruct_phase(R, s, wvl)
% R is reconstruction matrix
% s is slope vector (all x slopes then all y slopes)

OPD_m = R * s;  % OPD in meters
OPD = OPD_m * (2*pi) / wvl; % OPD in phase (radians)

% figure;
% plot(OPD_m .* 1e6, '-o');
% xlabel('OPD Index');
% ylabel('OPD (\mum)');

end

function [S, MSE] = Strehl_approximation(OPD1, OPD2)
% OPD1 & OPD2 are phase in radians

MSE = mean((OPD1 - OPD2).^2);
S = exp(-MSE);
end