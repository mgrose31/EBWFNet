clear; close all; clc;

%% load

filename_network = '.\Trained_Networks\LargeNetwork_SlopeMSE_CartesianMeshgrid_10msIE\depth4\trained_network_depth4_slopeMSE_cartesian_10msIE.mat';

% load the network
netH5 = load(filename_network, 'net');
net = netH5.net;

dd = dir('Formatted_Data\*.mat');

dataset_flag = 35;
for ii = 1:length(dd)
    if strfind(dd(ii).name, ['dataset', num2str(dataset_flag, '%02d')])
        idx_process = ii;
    end
end

folder_data = dd(idx_process).folder;
filename_data = dd(idx_process).name;
filepath_data = fullfile(folder_data, filename_data);

% load formatted file
tic; fprintf(1, "Loading formatted file " + filename_data + "... ");
dH5 = load(filepath_data, 'tvec_ev', 'positions', 'positions_baseline', ...
    'subImgs', 'xx', 'yy', 'dataset_flag', 'evFrame_ref', 'flirFrame_ref', ...
    'intertia_baseline', 'xc_ref_ev', 'yc_ref_ev', 'xc_ref_flir', ...
    'yc_ref_flir', 'tform_FLIR2ev');
fprintf(1, 'Done. '); toc;

positions_Frame = dH5.positions; % truth positions
predictions_AR1 = dH5.positions_baseline; % positions estimated by baseline inertia tracker

% convert time vector to seconds starting at 0
tvec_ev = (dH5.tvec_ev - dH5.tvec_ev(1)) ./ 1e6;

[numSubaps, numFrames] = size(dH5.subImgs, [4, 5]);

%% format data and make predictions

% X = dH5.subImgs(:, :, [1, 5], :, :); % depth of 1
% X = dH5.subImgs(:, :, [1, 2, 5, 6], :, :); % depth of 2
% X = dH5.subImgs(:, :, [1, 2, 3, 5, 6, 7], :, :); % depth of 3
X = dH5.subImgs; % full depth

xx = dH5.xx;
yy = dH5.yy;

% for mesh grids
xx = repmat(xx, [1, 1, 1, numSubaps, numFrames]);
yy = repmat(yy, [1, 1, 1, numSubaps, numFrames]);
X(:, :, end+1, :, :) = xx;
X(:, :, end+1, :, :) = yy;

% % for polar coordinates
% [theta, r] = cart2pol(xx, yy);
% % r = repmat(r, [1, 1, 1, numSubaps, numFrames]);
% % theta = repmat(theta, [1, 1, 1, numSubaps, numFrames]);
% X(:, :, end+1, :, :) = r;
% X(:, :, end+1, :, :) = theta;

X2 = dlarray(X(:,:,:,:), "SSCB"); %dlnetwork
clear('X');

%% run inference

minibatch_size = 1024;
num_inferences = size(X2, 4);
idx_inferences = 1:minibatch_size:num_inferences;
if idx_inferences(end) ~= num_inferences+1
    idx_inferences = [idx_inferences, num_inferences+1];
end

wb = waitbar(0, 'Making predictions... ');

predictions_CNN = zeros(2, num_inferences);

tstart = tic;
tend0_arr = zeros(1, length(idx_inferences)-1);
for ii = 1:length(idx_inferences)-1
    waitbar(idx_inferences(ii) / num_inferences, wb);

    X3 = gpuArray(X2(:,:,:,idx_inferences(ii):idx_inferences(ii+1)-1));

    % make the predictions
    % predictions_tmp = predict(net, X2(:,:,:,idx_inferences(ii):idx_inferences(ii+1)-1));
    predictions_tmp = predict(net, X3);

    % populate the array
    predictions_CNN(:, idx_inferences(ii):idx_inferences(ii+1)-1) = extractdata(predictions_tmp);

    tend0_arr(ii) = toc(tstart);
end
tend0 = toc(tstart);

predictions_CNN = reshape(predictions_CNN, [2, numSubaps, numFrames]);
fprintf(1, 'Done. '); toc;
close(wb);

clear('X3');

% [predictions_CNN(1,:,:), predictions_CNN(2,:,:)] = pol2cart(predictions_CNN(2,:,:), predictions_CNN(1,:,:));

%% quantify performance

% pix2rad = 5e-6; % 1 event pixel = 5 microrad of tilt
% 
% predictions_CNN_slope = predictions_CNN .* pix2rad;
% positions_FLIR_slope = positions_Frame .* pix2rad;
% positions_baseline_slope = predictions_AR1 .* pix2rad;

x_SE = squeeze((predictions_CNN(1, :, :) - positions_Frame(1, :, :)).^2);
y_SE = squeeze((predictions_CNN(2, :, :) - positions_Frame(2, :, :)).^2);
RSE_CNN = sqrt(x_SE + y_SE);

% x_SE_slope = squeeze((predictions_CNN_slope(1, :, :) - positions_FLIR_slope(1, :, :)).^2);
% y_SE_slope = squeeze((predictions_CNN_slope(2, :, :) - positions_FLIR_slope(2, :, :)).^2);
% RSE_CNN_slope = sqrt(x_SE_slope + y_SE_slope);

RSE_CNN_median = median(RSE_CNN(:));
RSE_CNN_mean = mean(RSE_CNN(:));

% RSE_CNN_slope_median = median(RSE_CNN_slope(:));
% RSE_CNN_slope_mean = mean(RSE_CNN_slope(:));

x_SE_baseline = squeeze((predictions_AR1(1, :, :) - positions_Frame(1, :, :)).^2);
y_SE_baseline = squeeze((predictions_AR1(2, :, :) - positions_Frame(2, :, :)).^2);
RSE_baseline = sqrt(x_SE_baseline + y_SE_baseline);

RSE_baseline_median = median(RSE_baseline(:));
RSE_baseline_mean = mean(RSE_baseline(:));

figure;
histogram(RSE_baseline, 0:0.1:12);
hold on;
histogram(RSE_CNN, 0:0.1:12);
xlim([0, 6]);
xlabel('root-squared-error');
ylabel('Counts');
legend('AR(1)', 'CNN (Ours)');
hold off;

%%

subap_plt = 20;

fig1 = figure('Position', [100, 100, 1000, 400]);
plot(tvec_ev, squeeze(positions_Frame(1, subap_plt, :)), 'LineWidth', 2);
hold on;
plot(tvec_ev, squeeze(predictions_AR1(1, subap_plt, :)), 'LineWidth', 2);
plot(tvec_ev, squeeze(predictions_CNN(1, subap_plt, :)), 'LineWidth', 2);
% xlim([tvec_ev(1), tvec_ev(end)]);
xlim([5, 6]);
xlabel('Time');
ylabel('x-axis spot position (event pixels)');
legend('Frame', 'AR(1)', 'CNN (Ours)', 'location', 'best');
set(gca, 'FontWeight', 'bold');
hold off;

% exportgraphics(fig1, 'Figure_xSlopes_dataset04_subaperture20.pdf');

fig2 = figure('Position', [100, 100, 1000, 400]);
plot(tvec_ev, squeeze(positions_Frame(2, subap_plt, :)), 'LineWidth', 2);
hold on;
plot(tvec_ev, squeeze(predictions_AR1(2, subap_plt, :)), 'LineWidth', 2);
plot(tvec_ev, squeeze(predictions_CNN(2, subap_plt, :)), 'LineWidth', 2);
% xlim([tvec_ev(1), tvec_ev(end)]);
xlim([5, 6]);
xlabel('Time');
ylabel('y-axis spot position (event pixels)');
legend('Frame', 'AR(1)', 'CNN (Ours)', 'location', 'best');
set(gca, 'FontWeight', 'bold');
hold off;

% exportgraphics(fig2, 'Figure_ySlopes_dataset04_subaperture20.pdf');

% figure('Position', [100, 100, 1000, 400]);
% subplot(1, 2, 1);
% binscatter(squeeze(positions_FLIR(1, subap_plt, :)), squeeze(predictions_CNN(1, subap_plt, :)), 100);
% hold on;
% h = plot([-20, 20], [-20, 20], 'k-', 'LineWidth', 1);
% xlim([-20, 20]);
% ylim([-20, 20]);
% legend(h, 'y=x', 'location', 'northwest');
% hold off;
% 
% subplot(1, 2, 2);
% binscatter(squeeze(positions_FLIR(2, subap_plt, :)), squeeze(predictions_CNN(2, subap_plt, :)), 100);
% hold on;
% h = plot([-20, 20], [-20, 20], 'k-', 'LineWidth', 1);
% xlim([-20, 20]);
% ylim([-20, 20]);
% legend(h, 'y=x', 'location', 'northwest');
% hold off;

[fx_truth, xx_truth] = ecdf(positions_Frame(1,:));
[fy_truth, yy_truth] = ecdf(positions_Frame(2,:));

[fx_baseline, xx_baseline] = ecdf(predictions_AR1(1,:));
[fy_baseline, yy_baseline] = ecdf(predictions_AR1(2,:));

[fx_CNN, xx_CNN] = ecdf(predictions_CNN(1,:));
[fy_CNN, yy_CNN] = ecdf(predictions_CNN(2,:));

S_CDF.fx_truth = fx_truth;
S_CDF.xx_truth = xx_truth;
S_CDF.fy_truth = fy_truth;
S_CDF.yy_truth = yy_truth;

S_CDF.fx_baseline = fx_baseline;
S_CDF.xx_baseline = xx_baseline;
S_CDF.fy_baseline = fy_baseline;
S_CDF.yy_baseline = yy_baseline;

S_CDF.fx_CNN = fx_CNN;
S_CDF.xx_CNN = xx_CNN;
S_CDF.fy_CNN = fy_CNN;
S_CDF.yy_CNN = yy_CNN;

% figure;
% plot(xx_truth, fx_truth, 'LineWidth', 2);
% hold on;
% plot(xx_baseline, fx_baseline, 'LineWidth', 2);
% plot(xx_CNN, fx_CNN, 'LineWidth', 2);
% xlabel('x-slope');
% ylabel('P[X \leq x]');
% legend('Truth', 'Baseline', 'CNN', 'location', 'southeast');
% hold off;
% 
% figure;
% plot(yy_truth, fy_truth, 'LineWidth', 2);
% hold on;
% plot(yy_baseline, fy_baseline, 'LineWidth', 2);
% plot(yy_CNN, fy_CNN, 'LineWidth', 2);
% xlabel('y-slope');
% ylabel('P[X \leq x]');
% legend('Truth', 'Baseline', 'CNN', 'location', 'southeast');
% hold off;

%%

% for subap_plt = 1:size(positions_truth, 2)
%     figure('Position', [100, 100, 1000, 400]);
%     plot(tvec_ev, squeeze(positions_truth(2, subap_plt, :)), 'LineWidth', 2);
%     hold on;
%     plot(tvec_ev, squeeze(predictions_NN(2, subap_plt, :)), 'LineWidth', 2);
%     plot(tvec_ev, squeeze(positions_baseline(2, subap_plt, :)), 'LineWidth', 2);
%     xlim([tvec_ev(1), tvec_ev(end)]);
%     xlabel('Time');
%     ylabel('y-axis centroid position');
%     title(['Sub-aperture ', num2str(subap_plt)]);
%     legend('Truth', 'CNN', 'Baseline', 'location', 'best');
%     hold off;
% end

%% saving

dataset_flag = dH5.dataset_flag;
evFrame_ref = dH5.evFrame_ref;
flirFrame_ref = dH5.flirFrame_ref;
inertia_baseline = dH5.intertia_baseline;
xc_ref_ev = dH5.xc_ref_ev;
yc_ref_ev = dH5.yc_ref_ev;
xc_ref_flir = dH5.xc_ref_flir;
yc_ref_flir = dH5.yc_ref_flir;
tform_FLIR2ev = dH5.tform_FLIR2ev;

save(['predictions_dataset', num2str(dataset_flag, '%02d'), '.mat'], ...
    'dataset_flag', 'evFrame_ref', 'flirFrame_ref', 'inertia_baseline', ...
    'xc_ref_ev', 'yc_ref_ev', 'xc_ref_flir', 'yc_ref_flir', 'tform_FLIR2ev', ...
    'positions_Frame', 'predictions_AR1', 'predictions_CNN', ...
    'RSE_CNN', 'RSE_baseline', 'RSE_CNN_mean', 'RSE_CNN_median', ...
    'RSE_baseline_mean', 'RSE_baseline_median', 'S_CDF', 'filename_network');

%%

% xx_plt = xx(1,:,1,1,1);
% yy_plt = yy(:,1,1,1,1);
% 
% subap_plt = 5;
% % frame_plt = 100;
% 
% v = VideoWriter('test_results_vid2', 'MPEG-4');
% v.Quality = 100;
% v.FrameRate = 5;
% 
% open(v);
% 
% fig1 = figure('Position', [100, 100, 1000, 400]);
% hold on;
% 
% % for frame_plt = 1:size(X, 5)
% for frame_plt = 1:500
% 
%     ax1 = subplot(1, 2, 1);
%     imagesc(xx_plt, yy_plt, X(:, :, 1, subap_plt, frame_plt));
%     axis image;
%     colormap(ax1, flipud(gray));
%     colorbar; clim([0, 10]);
%     hold on;
%     plot(positions_FLIR(1, subap_plt, frame_plt), ...
%         positions_FLIR(2, subap_plt, frame_plt), ...
%         'rx', 'LineWidth', 2, 'MarkerSize', 10);
%     plot(positions_baseline(1, subap_plt, frame_plt), ...
%         positions_baseline(2, subap_plt, frame_plt), ...
%         'bx', 'LineWidth', 2, 'MarkerSize', 10);
%     plot(predictions_CNN(1, subap_plt, frame_plt), ...
%         predictions_CNN(2, subap_plt, frame_plt), ...
%         'gx', 'LineWidth', 2, 'MarkerSize', 10);
%     title('Positive TORE Frame');
%     legend('Truth', 'Baseline', 'CNN', 'location', 'northeast');
%     hold off;
% 
%     ax2 = subplot(1, 2, 2);
%     imagesc(xx_plt, yy_plt, X(:, :, 3, subap_plt, frame_plt));
%     axis image;
%     colormap(ax2, flipud(gray));
%     colorbar; clim([0, 10]);
%     hold on;
%     plot(positions_FLIR(1, subap_plt, frame_plt), ...
%         positions_FLIR(2, subap_plt, frame_plt), ...
%         'rx', 'LineWidth', 2, 'MarkerSize', 10);
%     plot(positions_baseline(1, subap_plt, frame_plt), ...
%         positions_baseline(2, subap_plt, frame_plt), ...
%         'bx', 'LineWidth', 2, 'MarkerSize', 10);
%     plot(predictions_CNN(1, subap_plt, frame_plt), ...
%         predictions_CNN(2, subap_plt, frame_plt), ...
%         'gx', 'LineWidth', 2, 'MarkerSize', 10);
%     title('Negative TORE Frame');
%     legend('Truth', 'Baseline', 'CNN', 'location', 'northeast');
%     hold off;
% 
%     sgtitle(['Frame ', num2str(frame_plt)]);
% 
%     writeVideo(v, getframe(fig1));
%     clf(fig1);
% end
% 
% close(v);
% close(fig1);


%%
% 
% xx_plt = xx(1,:,1,1,1);
% yy_plt = yy(:,1,1,1,1);
% 
% subap_plt = 1;
% frame_plt = 100;
% 
% figure('Position', [100, 100, 1000, 400]);
% ax1 = subplot(1, 2, 1);
% imagesc(xx_plt, yy_plt, X(:, :, 1, subap_plt, frame_plt));
% axis image;
% colormap(ax1, flipud(gray));
% colorbar; clim([0, 10]);
% hold on;
% plot(positions_truth(1, subap_plt, frame_plt), ...
%     positions_truth(2, subap_plt, frame_plt), ...
%     'rx', 'LineWidth', 2, 'MarkerSize', 10);
% plot(positions_baseline(1, subap_plt, frame_plt), ...
%     positions_baseline(2, subap_plt, frame_plt), ...
%     'bx', 'LineWidth', 2, 'MarkerSize', 10);
% plot(predictions(1, subap_plt, frame_plt), ...
%     predictions(2, subap_plt, frame_plt), ...
%     'gx', 'LineWidth', 2, 'MarkerSize', 10);
% title('Positive TORE Frame');
% legend('Truth', 'Baseline', 'CNN', 'location', 'northeast');
% hold off;
% 
% ax2 = subplot(1, 2, 2);
% imagesc(xx_plt, yy_plt, X(:, :, 3, subap_plt, frame_plt));
% axis image;
% colormap(ax2, flipud(gray));
% colorbar; clim([0, 10]);
% hold on;
% plot(positions_truth(1, subap_plt, frame_plt), ...
%     positions_truth(2, subap_plt, frame_plt), ...
%     'rx', 'LineWidth', 2, 'MarkerSize', 10);
% plot(positions_baseline(1, subap_plt, frame_plt), ...
%     positions_baseline(2, subap_plt, frame_plt), ...
%     'bx', 'LineWidth', 2, 'MarkerSize', 10);
% plot(predictions(1, subap_plt, frame_plt), ...
%     predictions(2, subap_plt, frame_plt), ...
%     'gx', 'LineWidth', 2, 'MarkerSize', 10);
% title('Negative TORE Frame');
% legend('Truth', 'Baseline', 'CNN', 'location', 'northeast');
% hold off;