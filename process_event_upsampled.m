clear; close all; clc;
addpath('C:\Users\ISSL\Documents\ISSL Sync\Mitchell\Scripts');

dataset_flag = 4;

evFrames_FPS = 10000; % desired framerate of event data
FLIR_FPS = 100; % frame data collected at this rate

% load formatted data
load(['.\Formatted_Data\formatted_dataset', num2str(dataset_flag, '%02d'), '.mat'], ...
    'events', 'frameID', 'positions', 'positions_baseline', 'tvec_ev', ...
    'xc_ref_ev_rounded', 'yc_ref_ev_rounded', 'xc_ref_ev', 'yc_ref_ev');

% load prediction data
d_Predictions = load(['.\Trained_Networks\LargeNetwork_SlopeMSE_CartesianMeshgrid\depth4\predictions\predictions_dataset', num2str(dataset_flag, '%02d'), '.mat']);

positions_Frame = positions;
clear('positions');

inertia_baseline = d_Predictions.inertia_baseline;

%% set processing parameters

% upsample factor
n_upsample = evFrames_FPS / FLIR_FPS; % upsample factor

% height and width of event sensor
height = events.height;
width = events.width;

%%

tvec_ev2 = tvec_ev - tvec_ev(1);
tvec_ev3 = tvec_ev2 ./ 1e6;

idx_upsample = find(tvec_ev3 >= 13.70 & tvec_ev3 <= 14.00);

tvec_ev3_up = upsample(tvec_ev3(idx_upsample), n_upsample);
tvec_ev3_up = tvec_ev3_up(1:end-n_upsample+1);
idx_tmp = find(tvec_ev3_up);
tvec_ev4 = interp1(idx_tmp, tvec_ev3_up(idx_tmp), 1:length(tvec_ev3_up), "linear");

tvec_ev4 = tvec_ev4 .* 1e6;

events2 = events;
events2.ts = events2.ts - tvec_ev(1);

idx_keep = events2.ts >= 0;
events2.ts = events2.ts(idx_keep);
events2.x = events2.x(idx_keep);
events2.y = events2.y(idx_keep);
events2.p = events2.p(idx_keep);

evFrames_TORE_up = events2ToreFeature(events2.x, events2.y, events2.ts, events2.p, tvec_ev4, 4, [height, width]);

% % way to confirm things work; need to load original TORE frames
% all(evFrames_TORE_tmp(:,:,:,1:10:510) == evFrames_TORE(:,:,:,200:250), 'all')

%% set up images for inference

[xx, yy] = meshgrid(-25:25-1, -25:25-1); % define sub-aperture size

numFrames = size(evFrames_TORE_up, 4);
numSubaps = length(xc_ref_ev_rounded);

subImgs = NaN([size(xx), size(evFrames_TORE_up, 3) + 2, numSubaps, numFrames]);

for ii = 1:numFrames
    fprintf(1, ['Processing frame ', num2str(ii), ' of ', num2str(numFrames), '\n']);
    for jj = 1:length(xc_ref_ev_rounded)

        img_tmp = double(evFrames_TORE_up(yc_ref_ev_rounded(jj)-25:yc_ref_ev_rounded(jj)+25-1,...
            xc_ref_ev_rounded(jj)-25:xc_ref_ev_rounded(jj)+25-1,:,ii));

        img_tmp(:,:,end+1) = xx;
        img_tmp(:,:,end+1) = yy;
        subImgs(:,:,:,jj,ii) = img_tmp;

    end
end

frameID_up = 1:numFrames;

% clear('evFrames_TORE_up');

%% load and apply trained network

filename_network = './Trained_Networks/LargeNetwork_SlopeMSE_CartesianMeshgrid\depth4\trained_network_depth4_slopeMSE_cartesian.mat';

% load the network
netH5 = load(filename_network, 'net');
net = netH5.net;

X = dlarray(subImgs(:,:,:,:), "SSCB"); %dlnetwork

num_inferences = size(X, 4);
idx_inferences = [1:5000:num_inferences, num_inferences+1];

tic;
wb = waitbar(0, 'Making predictions... ');

predictions_CNN_up = zeros(2, num_inferences);
for ii = 1:length(idx_inferences)-1
    waitbar(idx_inferences(ii) / num_inferences, wb);
    X2 = gpuArray(X(:,:,:,idx_inferences(ii):idx_inferences(ii+1)-1));
    predictions_tmp = predict(net, X2);
    predictions_CNN_up(:, idx_inferences(ii):idx_inferences(ii+1)-1) = extractdata(predictions_tmp); % dlnetwork
end
predictions_CNN_up = reshape(predictions_CNN_up, [2, numSubaps, numFrames]);
fprintf(1, 'Done. '); toc;

close(wb);

clear('X');
clear('X2');

%% set plotting time vectors to start at time 0

tvec_ev_plot = (tvec_ev - tvec_ev(1)) ./ 1e6;
tvec_ev_plot_up = tvec_ev4 ./ 1e6;

%% paper AR1 autoregressive filter

m = inertia_baseline;

for ii = 1:length(xc_ref_ev)  % loop over each subaperture
    idx_events_subap = events2.y >= yc_ref_ev(ii) - 25 ...
        & events2.y <= yc_ref_ev(ii) + 24 ...
        & events2.x >= xc_ref_ev(ii) - 25 ...
        & events2.x <= xc_ref_ev(ii) + 24 ...
        & events2.p == 1;

    events_subap = events2;
    events_subap.x = events_subap.x(idx_events_subap);
    events_subap.y = events_subap.y(idx_events_subap);
    events_subap.p = events_subap.p(idx_events_subap);
    events_subap.ts = events_subap.ts(idx_events_subap);

    x_hat = NaN(size(events_subap.x));
    y_hat = NaN(size(events_subap.x));

    x_hat(1) = events_subap.x(1);
    y_hat(1) = events_subap.y(1);

    for jj = 2:length(events_subap.x)
        x_hat(jj) = m * x_hat(jj-1) + (1 - m) * events_subap.x(jj);
        y_hat(jj) = m * y_hat(jj-1) + (1 - m) * events_subap.y(jj);
    end

    % interpolate paper's tracks to FLIR times
    ts_hat = events_subap.ts;

    [C, ia, ~] = unique(ts_hat);

    x_hat_interp(ii, :) = interp1(C, x_hat(ia), tvec_ev4);
    y_hat_interp(ii, :) = interp1(C, y_hat(ia), tvec_ev4);

end

xc_ref_ev_rounded_mat = repmat(xc_ref_ev_rounded, [1, size(x_hat_interp, 2)]);
yc_ref_ev_rounded_mat = repmat(yc_ref_ev_rounded, [1, size(y_hat_interp, 2)]);

xc_hat_paper = x_hat_interp - xc_ref_ev_rounded_mat;
yc_hat_paper = y_hat_interp - yc_ref_ev_rounded_mat;

%% set up slope data for further processing

pix2rad = 5e-6; % 1 pixel is 5 microradians of tilt

predictions_AR1 = d_Predictions.positions_baseline;
predictions_CNN = d_Predictions.predictions_CNN; % 100 Hz predictions

sx_Frame = squeeze(positions_Frame(1,:,:)) .* pix2rad;
sy_Frame = squeeze(positions_Frame(2,:,:)) .* pix2rad;
sx_AR1 = squeeze(predictions_AR1(1,:,:)) .* pix2rad;
sy_AR1 = squeeze(predictions_AR1(2,:,:)) .* pix2rad;
sx_CNN = squeeze(predictions_CNN(1,:,:)) .* pix2rad;
sy_CNN = squeeze(predictions_CNN(2,:,:)) .* pix2rad;
sx_CNN_up = squeeze(predictions_CNN_up(1,:,:)) .* pix2rad;
sy_CNN_up = squeeze(predictions_CNN_up(2,:,:)) .* pix2rad;

sx_AR1_up = xc_hat_paper .* pix2rad;
sy_AR1_up = yc_hat_paper .* pix2rad;

%% quick plots

subap_plt = 20;

fig1 = figure('Position', [100, 100, 1000, 400]);
plot(tvec_ev_plot, sx_Frame(subap_plt,:) .* 1e6,'-o', 'LineWidth', 2);
hold on;
plot(tvec_ev_plot, sx_CNN(subap_plt,:) .* 1e6, '-o', 'LineWidth', 2);
plot(tvec_ev_plot_up, sx_CNN_up(subap_plt,:) .* 1e6, 'LineWidth', 2);
plot(tvec_ev_plot, sx_AR1(subap_plt,:) .* 1e6, '-o', 'LineWidth', 2);
plot(tvec_ev_plot_up, sx_AR1_up(subap_plt,:) .* 1e6, 'LineWidth', 2);
xlim([min(tvec_ev_plot_up), max(tvec_ev_plot_up)]);
% xlim([tvec_ev_up_plot(1), tvec_ev_up_plot(end)]);
% xlim([5, 6]);
xlabel('Time (s)');
ylabel('x-axis slopes (\murad)');
legend('FLIR 100 Hz', 'EBWFCNN 100 Hz', 'EBWFCNN 10000 Hz', 'AR(1) 100 Hz', 'AR(1) 10000 Hz', 'location', 'best');
hold off;

% exportgraphics(fig1, 'UpSampled_Tracks_xAxis.pdf');
% exportgraphics(fig1, 'UpSampled_Tracks_xAxis.png');

fig2 = figure('Position', [100, 100, 1000, 400]);
plot(tvec_ev_plot, sy_Frame(subap_plt,:) .* 1e6, '-o', 'LineWidth', 2);
hold on;
plot(tvec_ev_plot, sy_CNN(subap_plt,:) .* 1e6, '-o', 'LineWidth', 2);
plot(tvec_ev_plot_up, sy_CNN_up(subap_plt,:) .* 1e6, 'LineWidth', 2);
plot(tvec_ev_plot, sy_AR1(subap_plt,:) .* 1e6, '-o', 'LineWidth', 2);
plot(tvec_ev_plot_up, sy_AR1_up(subap_plt,:) .* 1e6, 'LineWidth', 2);
xlim([min(tvec_ev_plot_up), max(tvec_ev_plot_up)]);
% xlim([tvec_ev_up_plot(1), tvec_ev_up_plot(end)]);
% xlim([5, 6]);
xlabel('Time (s)');
ylabel('y-axis slopes (\murad)');
legend('FLIR 100 Hz', 'EBWFCNN 100 Hz', 'EBWFCNN 10000 Hz', 'AR(1) 100 Hz', 'AR(1) 10000 Hz', 'location', 'best');
hold off;

% exportgraphics(fig2, 'UpSampled_Tracks_yAxis.pdf');
% exportgraphics(fig2, 'UpSampled_Tracks_yAxis.png');



%% remove time-avg

sx_Frame_timeavg = mean(sx_Frame, 2);
sy_Frame_timeavg = mean(sy_Frame, 2);

sx_AR1_timeavg = mean(sx_AR1, 2);
sy_AR1_timeavg = mean(sy_AR1, 2);

sx_CNN_timeavg = mean(sx_CNN, 2);
sy_CNN_timeavg = mean(sy_CNN, 2);

% use the same time-avg as regularly sampled AR(1)
sx_AR1_up_timeavg = sx_AR1_timeavg;
sy_AR1_up_timeavg = sy_AR1_timeavg;

% use the same time-avg as regularly sampled CNN
sx_CNN_up_timeavg = sx_CNN_timeavg;
sy_CNN_up_timeavg = sy_CNN_timeavg;

sx_Frame = sx_Frame - repmat(sx_Frame_timeavg, [1, size(sx_Frame, 2)]);
sy_Frame = sy_Frame - repmat(sy_Frame_timeavg, [1, size(sy_Frame, 2)]);

sx_AR1 = sx_AR1 - repmat(sx_AR1_timeavg, [1, size(sx_AR1, 2)]);
sy_AR1 = sy_AR1 - repmat(sy_AR1_timeavg, [1, size(sy_AR1, 2)]);

sx_AR1_up = sx_AR1_up - repmat(sx_AR1_up_timeavg, [1, size(sx_AR1_up, 2)]);
sy_AR1_up = sy_AR1_up - repmat(sy_AR1_up_timeavg, [1, size(sy_AR1_up, 2)]);

sx_CNN = sx_CNN - repmat(sx_CNN_timeavg, [1, size(sx_CNN, 2)]);
sy_CNN = sy_CNN - repmat(sy_CNN_timeavg, [1, size(sy_CNN, 2)]);

sx_CNN_up = sx_CNN_up - repmat(sx_CNN_up_timeavg, [1, size(sx_CNN_up, 2)]);
sy_CNN_up = sy_CNN_up - repmat(sy_CNN_up_timeavg, [1, size(sy_CNN_up, 2)]);

%% remove full-aperture tilt

sx_Frame_fullApTilt = mean(sx_Frame, 1);
sy_Frame_fullApTilt = mean(sy_Frame, 1);

sx_AR1_fullApTilt = mean(sx_AR1, 1);
sy_AR1_fullApTilt = mean(sy_AR1, 1);

sx_AR1_up_fullApTilt = mean(sx_AR1_up, 1);
sy_AR1_up_fullApTilt = mean(sy_AR1_up, 1);

sx_CNN_fullApTilt = mean(sx_CNN, 1);
sy_CNN_fullApTilt = mean(sy_CNN, 1);

sx_CNN_up_fullApTilt = mean(sx_CNN_up, 1);
sy_CNN_up_fullApTilt = mean(sy_CNN_up, 1);

sx_Frame = sx_Frame - repmat(sx_Frame_fullApTilt, [size(sx_Frame, 1), 1]);
sy_Frame = sy_Frame - repmat(sy_Frame_fullApTilt, [size(sy_Frame, 1), 1]);

sx_AR1 = sx_AR1 - repmat(sx_AR1_fullApTilt, [size(sx_AR1, 1), 1]);
sy_AR1 = sy_AR1 - repmat(sy_AR1_fullApTilt, [size(sy_AR1, 1), 1]);

sx_AR1_up = sx_AR1_up - repmat(sx_AR1_up_fullApTilt, [size(sx_AR1_up, 1), 1]);
sy_AR1_up = sy_AR1_up - repmat(sy_AR1_up_fullApTilt, [size(sy_AR1_up, 1), 1]);

sx_CNN = sx_CNN - repmat(sx_CNN_fullApTilt, [size(sx_CNN, 1), 1]);
sy_CNN = sy_CNN - repmat(sy_CNN_fullApTilt, [size(sy_CNN, 1), 1]);

sx_CNN_up = sx_CNN_up - repmat(sx_CNN_up_fullApTilt, [size(sx_CNN_up, 1), 1]);
sy_CNN_up = sy_CNN_up - repmat(sy_CNN_up_fullApTilt, [size(sy_CNN_up, 1), 1]);

%% plot

subap_plt = 7;

% figure('Position', [100, 100, 1000, 200]);
% plot(tvec_ev_plot, sx_Frame(subap_plt,:) .* 1e6, 'k-*', 'LineWidth', 2);
% hold on;
% plot(tvec_ev_plot, sx_AR1(subap_plt,:) .* 1e6, '-*', 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
% plot(tvec_ev_plot, sx_CNN(subap_plt,:) .* 1e6, '-*', 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
% plot(tvec_ev_plot_up, sx_CNN_up(subap_plt,:) .* 1e6, '-', 'LineWidth', 2, 'Color', [0.9290 0.6940 0.1250]);
% % xlim([11, 14]);
% xlim([13.80, 13.95]);
% ylim([-25, 25]);
% xlabel('Time (s)');
% ylabel('x-axis slope (\murad)');
% legend('Frame', 'AR(1)', 'EBWFNet (100 Hz)', 'EBWFNet (1000 Hz)', 'location', 'best');
% hold off;
% 
% figure('Position', [100, 100, 1000, 200]);
% plot(tvec_ev_plot, sy_Frame(subap_plt,:) .* 1e6, 'k-', 'LineWidth', 2);
% hold on;
% plot(tvec_ev_plot, sy_AR1(subap_plt,:) .* 1e6, '-*', 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
% plot(tvec_ev_plot, sy_CNN(subap_plt,:) .* 1e6, '-*', 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
% plot(tvec_ev_plot_up, sy_CNN_up(subap_plt,:) .* 1e6, '-', 'LineWidth', 2, 'Color', [0.9290 0.6940 0.1250]);
% % xlim([11, 14]);
% xlim([13.80, 13.95]);
% ylim([-25, 25]);
% xlabel('Time (s)');
% ylabel('y-axis slope (\murad)');
% legend('Frame', 'AR(1)', 'EBWFNet (100 Hz)', 'EBWFNet (1000 Hz)', 'location', 'best');
% hold off;


fig_xSlopes_upsampled = figure('Position', [100, 100, 1000, 200]);
plot(tvec_ev_plot, sx_Frame(subap_plt,:) .* 1e6, 'k-*', 'LineWidth', 2);
hold on;
plot(tvec_ev_plot_up, sx_CNN_up(subap_plt,:) .* 1e6, '-', 'LineWidth', 2, 'Color', [0.9290 0.6940 0.1250]);
xlim([13.80, 13.90]);
% ylim([-25, 25]);
% ylim([-10, 10]);
xlabel('Time (s)');
ylabel('x-axis slope (\murad)');
legend('Frame (100 Hz)', 'EBWFNet (10000 Hz)', 'location', 'best');
set(gca, 'FontWeight', 'bold');
hold off;

fig_ySlopes_upsampled = figure('Position', [100, 100, 1000, 200]);
plot(tvec_ev_plot, sy_Frame(subap_plt,:) .* 1e6, 'k-o', 'LineWidth', 2);
hold on;
plot(tvec_ev_plot_up, sy_CNN_up(subap_plt,:) .* 1e6, '-', 'LineWidth', 2, 'Color', [0.9290 0.6940 0.1250]);
xlim([13.80, 13.90]);
% ylim([-25, 25]);
% ylim([-10, 10]);
xlabel('Time (s)');
ylabel('y-axis slope (\murad)');
legend('Frame (100 Hz)', 'EBWFNet (10000 Hz)', 'location', 'best');
set(gca, 'FontWeight', 'bold');
hold off;

% exportgraphics(fig_xSlopes_upsampled, 'Figure_xSlopes_upsampled.pdf');
% exportgraphics(fig_ySlopes_upsampled, 'Figure_ySlopes_upsampled.pdf');

%%

fig_Slopes_Upsampled_3D_View1 = figure;
scatter3(sx_Frame(subap_plt,:) .* 1e6, sy_Frame(subap_plt,:) .* 1e6, tvec_ev_plot, 'ko', 'MarkerFaceColor', 'k');
hold on;
plot3(sx_AR1_up(subap_plt,:) .* 1e6, sy_AR1_up(subap_plt,:) .* 1e6, tvec_ev_plot_up, '-', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 3);
plot3(sx_CNN_up(subap_plt,:) .* 1e6, sy_CNN_up(subap_plt,:) .* 1e6, tvec_ev_plot_up, '-', 'Color', [0.8500, 0.3250, 0.0980],  'LineWidth', 3);
xlim([-25, 25]);
ylim([-25, 25]);
zlim([13.75, 13.95]);
xticks(-30:10:30);
yticks(-30:10:30);
zticks(13.75:0.05:13.95);
xlabel('x-slope (\murad)');
ylabel('y-slope (\murad)');
zlabel('Time (s)');
legend('Frame', 'AR(1)', 'EBWFNet', 'location', 'best');
grid on;
% box on;
set(gca, 'FontWeight', 'bold', 'FontSize', 10);
hold off;


fig_Slopes_Upsampled_3D_View2 = figure;
scatter3(sx_Frame(subap_plt,:) .* 1e6, sy_Frame(subap_plt,:) .* 1e6, tvec_ev_plot, 'ko', 'MarkerFaceColor', 'k');
hold on;
plot3(sx_AR1_up(subap_plt,:) .* 1e6, sy_AR1_up(subap_plt,:) .* 1e6, tvec_ev_plot_up, '-', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 3);
plot3(sx_CNN_up(subap_plt,:) .* 1e6, sy_CNN_up(subap_plt,:) .* 1e6, tvec_ev_plot_up, '-', 'Color', [0.8500, 0.3250, 0.0980],  'LineWidth', 3);
xlim([-25, 25]);
ylim([-25, 25]);
zlim([13.75, 13.95]);
view(135, 10);
xticks(-30:10:30);
yticks(-30:10:30);
zticks(13.75:0.05:13.95);
xlabel('x-slope (\murad)');
ylabel('y-slope (\murad)');
zlabel('Time (s)');
legend('Frame', 'AR(1)', 'EBWFNet', 'location', 'best');
grid on;
% box on;
set(gca, 'FontWeight', 'bold', 'FontSize', 12);
hold off;

% exportgraphics(fig_Slopes_Upsampled_3D_View1, 'Figure_Slopes_Upsampled_3D_View1.pdf');
% exportgraphics(fig_Slopes_Upsampled_3D_View2, 'Figure_Slopes_Upsampled_3D_View2.pdf');

%% make a video

subap_plt = 7;
frame_plt = 100;

cmin = min(subImgs(:,:,1,:,1:10), [], 'all');
cmax = max(subImgs(:,:,1,:,1:10), [], 'all');

% figure;
% imagesc(xx(1,:), yy(:,1), subImgs(:, :, 1, subap_plt, frame_plt));
% axis image xy off;
% colormap(flipud(gray));
% colorbar; clim([cmin, cmax]);
% hold on;
% h1 = plot(predictions_CNN_up(1, subap_plt, frame_plt), predictions_CNN_up(2, subap_plt, frame_plt), 'gx', 'MarkerSize', 10, 'LineWidth', 3);
% h2 = plot(xc_hat_paper(subap_plt, frame_plt), yc_hat_paper(subap_plt, frame_plt), 'rx', 'MarkerSize', 10, 'LineWidth', 3);
% legend([h1, h2], {'AR(1)', 'EBWFNet'}, 'location', 'northeast');
% set(gca, 'FontWeight', 'Bold');
% hold off;

idx_tvec_ev_plot = find(tvec_ev_plot_up >= 13.75 & tvec_ev_plot_up <= 13.95);

flag_mkVid = 1;
if flag_mkVid
    v = VideoWriter('Tracks_Upsampled_test2', 'MPEG-4');
    v.Quality = 100;
    v.FrameRate = 50;
    open(v);
end

fig = figure;
hold on;
% for ii = 1:1:300
% for ii = 1:1:length(tvec_ev_plot_up)
for ii = idx_tvec_ev_plot
    frame_plt = ii;

    % imagesc(xx(1,:), yy(:,1), subImgs(:, :, 1, subap_plt, frame_plt));  % positive TORE frame depth K = 1
    imagesc(xx(1,:), yy(:,1), subImgs(:, :, 5, subap_plt, frame_plt));  % negative TORE frame depth K = 1
    axis image xy off;
    colormap(flipud(gray));
    colorbar; clim([cmin, cmax]);
    hold on;
    h1 = plot(xc_hat_paper(subap_plt, frame_plt), yc_hat_paper(subap_plt, frame_plt), 'gx', 'MarkerSize', 10, 'LineWidth', 3);
    h2 = plot(predictions_CNN_up(1, subap_plt, frame_plt), predictions_CNN_up(2, subap_plt, frame_plt), 'rx', 'MarkerSize', 10, 'LineWidth', 3);
    % title(['Frame ', num2str(frame_plt, '%04d')]);
    title(['Time: ', num2str(tvec_ev_plot_up(ii), '%.4f')]);
    legend([h1, h2], {'AR(1)', 'EBWFNet'}, 'location', 'northeast');
    set(gca, 'FontWeight', 'Bold');
    hold off;

    drawnow;

    if flag_mkVid
        writeVideo(v, getframe(fig));
    end
    clf(fig);
end
close(fig);
if flag_mkVid
    close(v);
end
