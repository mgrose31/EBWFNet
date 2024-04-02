clear; close all; clc;
% Code to format data for CNN training. Autoregressive filter from prior
% paper (AR1) is computed on event data. Transformation of frame (FLIR)
% tracks to event camera space is computed. Data is formatted and saved for
% training the CNN.
% Author: Mitchell Grose, University of Dayton, 2023

%% script hyperparameters

% dataset to process
idx_process = 35;

% subaperture sizes
subap_sz_ev = 54; % Prophesee VGA sensor
subap_sz_flir = 200; % FLIR frame sensor

%% load data

dd_ref_flir = dir(fullfile(pwd, './spot_ref_centers_20230730/*flir*.mat'));
dd_ref_ev = dir(fullfile(pwd, './spot_ref_centers_20230730/*ev*.mat'));
dd_tracks_flir = dir(fullfile(pwd, './Tracks_FLIR/*flir*.mat'));
dd_evFrames_TORE = dir(fullfile('D:/TORE_Volumes_20230730_02p5msIE/*.mat'));

str_find = ['dataset', num2str(idx_process, '%02d')];

for ii=1:length(dd_evFrames_TORE)
    if contains(dd_evFrames_TORE(ii).name, str_find)
        idx_evFrames_TORE = ii;
    end
end

for ii=1:length(dd_ref_ev)
    if contains(dd_ref_ev(ii).name, str_find)
        idx_ref_ev = ii;
    end
end

for ii=1:length(dd_ref_flir)
    if contains(dd_ref_flir(ii).name, str_find)
        idx_ref_flir = ii;
    end
end

for ii=1:length(dd_tracks_flir)
    if contains(dd_tracks_flir(ii).name, str_find)
        idx_tracks_flir = ii;
    end
end

%% load event frames

tic; fprintf(1, 'Loading event frames... ');
load(fullfile(dd_evFrames_TORE(idx_evFrames_TORE).folder,...
    dd_evFrames_TORE(idx_evFrames_TORE).name));
fprintf(1, 'Done. '); toc;

%% set frames to keep
% Not all frames are kept for a variety of reasons. Sometimes a car drives
% by and interrupts the data. Other times the subapertures run off the edge
% of the frame camera. First 200 frames are always suspect due to
% misordering of frames by Spinview, so we always exclude first 200.

numFrames = size(evFrames_TORE, 4);

if ismember(idx_process, [1:3, 5:15, 17, 24, 25, 26, 27, 32, 34])
    idx_keep_frames = [200:900, 1100:numFrames];
elseif idx_process == 4
    idx_keep_frames = [200:900, 1100:1850]; % exclude specific frames
elseif idx_process == 16
    idx_keep_frames = [200:900, 1100:2100]; % exclude frames after 2200
elseif idx_process == 19
    idx_keep_frames = [700:900, 1100:numFrames];
elseif idx_process == 22
    idx_keep_frames = [200:900, 1300:numFrames]; % exclude 1100-1300 because signal disappears momentarily
elseif idx_process == 28
    idx_keep_frames = [200:900, 1100:2400];
elseif idx_process == 29
    idx_keep_frames = 200:900;
elseif idx_process == 30
    idx_keep_frames = [200:900, 1500:numFrames]; % exclude 1100:1499 because car drives by
elseif idx_process == 31
    idx_keep_frames = [500:900, 1100:1700]; % exclude 200:500 because spots between subapertures
elseif idx_process == 33
    idx_keep_frames = [200:900, 1100:1800]; % spots move starting after 1800
elseif idx_process == 35
    idx_keep_frames = [200:900, 1100:2450]; % ignore frames at end because car drives by
end

%% load event reference positions

% load file with reference subaperture positions
fpath_ev_subapRefs = fullfile(dd_ref_ev(idx_ref_ev).folder,...
    dd_ref_ev(idx_ref_ev).name);

% Event sensor
if isfile(fpath_ev_subapRefs)
    load(fpath_ev_subapRefs);
else
    error('File does not exist!');
end

% rectangles to draw for each subaperture position
r_pos_ev = [xc_ref_ev-subap_sz_ev/2, yc_ref_ev-subap_sz_ev/2,...
    ones(size(yc_ref_ev)).*subap_sz_ev, ones(size(yc_ref_ev)).*subap_sz_ev];

load('positions_hot_pixels.mat');
% evFrame_ref(positions_hot_pixel) = 0;

idx_hot = zeros(size(events.ts));

for ii=1:length(positions_hot_pixel)
    idx_hot_tmp = events.x == positions_hot_pixel(ii,1) & events.y == positions_hot_pixel(ii,2);
    idx_hot = idx_hot | idx_hot_tmp;
end

events.x = events.x(~idx_hot);
events.y = events.y(~idx_hot);
events.p = events.p(~idx_hot);
events.ts = events.ts(~idx_hot);

evFrame_ref = accumarray([events.y, events.x], 1, [events.height, events.width], @sum, 0);

% plot reference image with subaperture positions
figure;
imagesc(evFrame_ref);
hold on;
plot(xc_ref_ev, yc_ref_ev, 'r+', 'MarkerSize', 10);
for ii=1:length(xc_ref_ev)
    rectangle('Position', r_pos_ev(ii,:), 'EdgeColor', 'g', 'LineWidth', 3);
end
axis image xy;
colorbar;
colormap('gray');
title('Compute Subaperture Reference Positions');
hold off;
drawnow;

%% load flir reference positions and tracks

fpath_flir_subapRefs = fullfile(dd_ref_flir(idx_ref_flir).folder,...
    dd_ref_flir(idx_ref_flir).name);
fpath_flir_tracks = fullfile(dd_tracks_flir(idx_tracks_flir).folder,...
    dd_tracks_flir(idx_tracks_flir).name);

if isfile(fpath_flir_subapRefs)
    load(fpath_flir_subapRefs)
else
    error('Path does not exist!');
end

if isfile(fpath_flir_tracks)
    load(fpath_flir_tracks);
else
    error('Path does not exist!');
end

% define rectangles for subaperture positions
r_pos_flir = [xc_ref_flir-subap_sz_flir/2, yc_ref_flir-subap_sz_flir/2,...
    ones(size(yc_ref_flir)).*subap_sz_flir, ones(size(yc_ref_flir)).*subap_sz_flir];

% plot reference image with subaperture positions
figure;
imagesc(flirFrame_ref);
hold on;
plot(xc_ref_flir, yc_ref_flir, 'r+', 'MarkerSize', 10);
for ii=1:length(xc_ref_flir)
    rectangle('Position', r_pos_flir(ii,:), 'EdgeColor', 'g', 'LineWidth', 3);
end
axis image xy;
colorbar;
colormap('gray');
title('Compute Subaperture Reference Positions');
hold off;
drawnow;

%% paper inertia tracker (autoregressive filter, AR(1))

% inertia of prior estimate
% m_arr = [0.01:0.01:0.98, 0.981:0.001:0.999];
m_arr = 0.98:0.001:0.999;

% loop over inertias
for kk = 1:length(m_arr)

    disp(kk); % display current index

    m = m_arr(kk); % set current inertia

    % loop over subapertures
    for ii = 1:size(xc_imgs_flir, 1)
        % find events in current subaperture
        idx_events_subap = events.y >= yc_ref_ev(ii) - 25 ...
            & events.y <= yc_ref_ev(ii) + 24 ...
            & events.x >= xc_ref_ev(ii) - 25 ...
            & events.x <= xc_ref_ev(ii) + 24 ...
            & events.p == 1;

        % create event stream for current subaperture
        events_subap = events;
        events_subap.x = events_subap.x(idx_events_subap);
        events_subap.y = events_subap.y(idx_events_subap);
        events_subap.p = events_subap.p(idx_events_subap);
        events_subap.ts = events_subap.ts(idx_events_subap);

        % pre-allocate x-axis and y-axis estimates
        x_hat = NaN(size(events_subap.x));
        y_hat = NaN(size(events_subap.x));

        % first estimate is the first event
        x_hat(1) = events_subap.x(1);
        y_hat(1) = events_subap.y(1);

        % loop over events and compute estimate
        for jj = 2:length(events_subap.x)
            x_hat(jj) = m * x_hat(jj-1) + (1 - m) * events_subap.x(jj);
            y_hat(jj) = m * y_hat(jj-1) + (1 - m) * events_subap.y(jj);
        end

        % get times of AR(1) estimates
        ts_hat = events_subap.ts;

        % must only consider unique times
        [C, ia, ~] = unique(ts_hat);

        % interpolate AR(1) estimates to frame times
        x_hat_interp(ii, :, kk) = interp1(C, x_hat(ia), tvec_ev);
        y_hat_interp(ii, :, kk) = interp1(C, y_hat(ia), tvec_ev);

    end
end

%% fit transformation from frame space to event space (first time)

% subaperture reference positions
x1 = xc_ref_flir;
y1 = yc_ref_flir;
x2 = xc_ref_ev;
y2 = yc_ref_ev;

% transform frame subaperture reference positions to event space
tform_FLIR2ev = fitgeotform2d([x1, y1], [x2, y2], 'projective');

% transform frame spot tracks to event space
[xc_imgs_flirT, yc_imgs_flirT] = transformPointsForward(tform_FLIR2ev, xc_imgs_flir, yc_imgs_flir);

% loop over AR(1) inertia values
for ii = 1:length(m_arr)

    % get the data at frame times
    x_hat_interp_use = x_hat_interp(:, :, ii);
    y_hat_interp_use = y_hat_interp(:, :, ii);

    % compute root-squared-error between AR(1) and transformed frame tracks
    x_SE_paper = (xc_imgs_flirT - x_hat_interp_use).^2;
    y_SE_paper = (yc_imgs_flirT - y_hat_interp_use).^2;
    RSE_paper = sqrt(x_SE_paper + y_SE_paper);

    % compute the median and mean error between AR(1) and frame estimates
    RSE_paper_median(ii) = median(RSE_paper, 'all', 'omitnan');
    RSE_paper_mean(ii) = mean(RSE_paper, 'all', 'omitnan');

end

% plot the median and mean error as a function of inertia
% looking for the minimum of the curve
figure;
plot(m_arr, RSE_paper_median, 'LineWidth', 2);
hold on;
plot(m_arr, RSE_paper_mean, 'LineWidth', 2);
xlabel('intertia');
ylabel('root-squared-error');
title('First Iteration');
legend('median', 'mean');
hold off;
drawnow;

% optimal inertia is the one that minimizes the mean error
% [~, I_m_optimal] = min(RSE_paper_median);
[~, I_m_optimal] = min(RSE_paper_mean);
m_optimal = m_arr(I_m_optimal);

% temporary optimal AR(1) estimate is the one that mimimizes mean error
xc_hat_paper_tmp = x_hat_interp(:, :, I_m_optimal);
yc_hat_paper_tmp = y_hat_interp(:, :, I_m_optimal);

%% fit transformation from frame space to event space (second time)

% use mean position of tracks in each sensor and axis
% only use valid frames
xc_imgs_flir_mean = mean(xc_imgs_flir(:, idx_keep_frames), 2, 'omitnan');
yc_imgs_flir_mean = mean(yc_imgs_flir(:, idx_keep_frames), 2, 'omitnan');
xc_hat_paper_mean = mean(xc_hat_paper_tmp(:, idx_keep_frames), 2, 'omitnan');
yc_hat_paper_mean = mean(yc_hat_paper_tmp(:, idx_keep_frames), 2, 'omitnan');

% set variables to compute transformation
x1 = xc_imgs_flir_mean;
y1 = yc_imgs_flir_mean;
x2 = xc_hat_paper_mean;
y2 = yc_hat_paper_mean;

% compute transform from frame space to event space
tform_FLIR2ev = fitgeotform2d([x1, y1], [x2, y2], 'projective');

% transform frame average positions to event space
[x1T, y1T] = transformPointsForward(tform_FLIR2ev, x1, y1);

% transform frame tracks to event space
[xc_imgs_flirT, yc_imgs_flirT] = transformPointsForward(tform_FLIR2ev, xc_imgs_flir, yc_imgs_flir);

% plot event frame, average event spots, and average frame spots
% looking for good overlap between average event spots and frame spots
figure;
imagesc(evFrame_ref);
axis image;
colormap(gray); colorbar;
hold on;
% plot(xc_ref_ev, yc_ref_ev, 'gx', 'LineWidth', 2);
plot(x2, y2, 'gx', 'LineWidth', 2);
plot(x1T, y1T, 'rx', 'LineWidth', 2);
legend('Event Spot Centers', 'Transformed FLIR Spot Centers');
hold off;

% loop over inertias
for ii = 1:length(m_arr)
    % get the data at frame times
    x_hat_interp_use = x_hat_interp(:, :, ii);
    y_hat_interp_use = y_hat_interp(:, :, ii);

    % compute root-squared-error between AR(1) and transformed frame tracks
    x_SE_paper = (xc_imgs_flirT - x_hat_interp_use).^2;
    y_SE_paper = (yc_imgs_flirT - y_hat_interp_use).^2;
    RSE_paper = sqrt(x_SE_paper + y_SE_paper);

    % get the median and mean error
    RSE_paper_median(ii) = median(RSE_paper, 'all', 'omitnan');
    RSE_paper_mean(ii) = mean(RSE_paper, 'all', 'omitnan');

end

% plot the median and mean error vs. inertia
figure;
plot(m_arr, RSE_paper_median, 'LineWidth', 2);
hold on;
plot(m_arr, RSE_paper_mean, 'LineWidth', 2);
xlabel('intertia');
ylabel('root-squared-error');
title('Second Iteration');
legend('median', 'mean');
hold off;
drawnow;

% optimal inertia is the one that minimizes the mean error
% [~, I_m_optimal] = min(RSE_paper_median);
[~, I_m_optimal] = min(RSE_paper_mean);
m_optimal = m_arr(I_m_optimal);

% optimal AR(1) estimate is the one that mimimizes mean error
xc_hat_paper = x_hat_interp(:, :, I_m_optimal);
yc_hat_paper = y_hat_interp(:, :, I_m_optimal);

%% confirm transformation

% compute the difference in mean spot position between frame tracks
% transformed into event space and the AR(1) tracks
xdiff = mean(xc_imgs_flirT(:, idx_keep_frames), 2, 'omitnan') - mean(xc_hat_paper(:, idx_keep_frames), 2, 'omitnan');
ydiff = mean(yc_imgs_flirT(:, idx_keep_frames), 2, 'omitnan') - mean(yc_hat_paper(:, idx_keep_frames), 2, 'omitnan');

% plot the differences as a quiver
% a random pattern (with small magnitude) indicates good transformation
figure;
quiver(xc_ref_ev, yc_ref_ev, xdiff, ydiff);

numSubaps = size(xc_hat_paper, 1);

% show the overlaid tracks (frame tracks transformed to event space, AR(1))
fig1 = figure('Position', [100, 100, 1000, 400]);
hold on;
for ii = 1:numSubaps
    plot(tvec_ev, xc_imgs_flirT(ii,:), 'LineWidth', 2);
    hold on;
    plot(tvec_ev, xc_hat_paper(ii,:), 'LineWidth', 2);
    xlim([tvec_ev(1), tvec_ev(end)]);
    xlabel('Time');
    ylabel('x-slope');
    title(['Subaperture ', num2str(ii, '%02d')]);
    drawnow;

    pause(1);
    clf(fig1);
end
close(fig1);

%% compute spot positions relative to corresponding subaperture

% round the event reference positions
xc_ref_ev_rounded = round(xc_ref_ev);
yc_ref_ev_rounded = round(yc_ref_ev);

% convert rounded event reference positions to a matrix
xc_ref_ev_rounded_mat = repmat(xc_ref_ev_rounded, [1, size(evFrames_TORE, 4)]);
yc_ref_ev_rounded_mat = repmat(yc_ref_ev_rounded, [1, size(evFrames_TORE, 4)]);

% center the frame tracks (in event space) to each subaperture
% this makes the spot displacement relative to each subaperture, not the
% entire frame
xc_imgs_flirT_new = xc_imgs_flirT - xc_ref_ev_rounded_mat;
yc_imgs_flirT_new = yc_imgs_flirT - yc_ref_ev_rounded_mat;

% center the event tracks to each subaperture
% this makes the spot displacement relative to each subaperture, not the
% entire frame
xc_hat_paper_new = xc_hat_paper - xc_ref_ev_rounded_mat;
yc_hat_paper_new = yc_hat_paper - yc_ref_ev_rounded_mat;

% populate variables with positions for further processing
% this is for AR(1) filter
positions_baseline(1,:,:) = xc_hat_paper_new;
positions_baseline(2,:,:) = yc_hat_paper_new; 
intertia_baseline = m_optimal;

%% create subaperture TORE volumes from entire event sensor

[xx, yy] = meshgrid(-25:25-1, -25:25-1); % subaperture pixel map meshgrid

% number of frames to process
NumFrames = size(evFrames_TORE, 4);

% preallocate subaperture images (TORE volumes)
subImgs = NaN([size(xx), size(evFrames_TORE, 3), length(xc_ref_ev), NumFrames]);

% frame positions in event space (proxy truth for CNN training)
positions = NaN([2, length(xc_ref_ev), NumFrames]);

% loop over frames and subapertures
cnt = 0;
for ii = 1:NumFrames
    fprintf(1, ['Processing frame ', num2str(ii), ' of ', num2str(NumFrames), '\n']);
    for jj = 1:length(xc_ref_ev)

        % get the TORE volume for current frame and subaperture
        img_tmp = double(evFrames_TORE(yc_ref_ev_rounded(jj)-25:yc_ref_ev_rounded(jj)+25-1,...
            xc_ref_ev_rounded(jj)-25:xc_ref_ev_rounded(jj)+25-1,:,ii));

        % populate subaperture image
        subImgs(:,:,:,jj,ii) = img_tmp;

        % populate proxy truth spot displacement
        positions(1,jj,ii) = xc_imgs_flirT_new(jj,ii);
        positions(2,jj,ii) = yc_imgs_flirT_new(jj,ii);

    end
end

% define vector of frames
frameID = 1:NumFrames;

%% find any frames with NaN

frame_hasnan = squeeze(any(isnan(positions), [1, 2]));
idx_frame_isnan = find(frame_hasnan);

if ismember(idx_frame_isnan, idx_keep_frames)
    idx_tmp = find(idx_frame_isnan == idx_keep_frames);
    idx_keep_frames(idx_tmp) = [];
end

%% keep only good frames

frameID = frameID(idx_keep_frames);
positions = positions(:, :, idx_keep_frames);
subImgs = subImgs(:, :, :, :, idx_keep_frames);
tvec_ev = tvec_ev(idx_keep_frames);
positions_baseline = positions_baseline(:, :, idx_keep_frames);

%% save the data

% save the data as a single formatted file
tic; fprintf(1, 'Saving formatted file... ');
save(['formatted_dataset', num2str(dataset_flag, '%02d'), '.mat'], ...
    'frameID', 'positions', 'subImgs', 'tvec_ev', 'xx', 'yy', ...
    'tform_FLIR2ev', 'xc_ref_ev', 'yc_ref_ev', 'xc_ref_flir', ...
    'yc_ref_flir', 'xc_ref_ev_rounded', 'yc_ref_ev_rounded', ...
    'evFrame_ref', 'flirFrame_ref', 'dataset_flag', ...
    'positions_baseline', 'intertia_baseline', 'events', ...
    'xdiff', 'ydiff', 'idx_frame_isnan', '-v7.3');
fprintf(1, 'Saving individual files... ');


% save individual files
xx_tmp = repmat(xx, [1, 1, 1, size(subImgs, [4, 5])]);
yy_tmp = repmat(yy, [1, 1, 1, size(subImgs, [4, 5])]);

% append subaperture TORE volumes with xx and yy pixel maps (meshgrid)
subImgs(:,:,end+1,:,:) = xx_tmp;
subImgs(:,:,end+1,:,:) = yy_tmp;

% make save folder if it does not exist
% folder_save = fullfile(pwd, ['dataset', num2str(dataset_flag, '%02d')]);
folder_save = fullfile('C:\Users\ISSL\Documents\Formatted_Data_CNN_2p5msIE\', ['dataset', num2str(dataset_flag, '%02d')]);
if ~isfolder(folder_save)
    mkdir(folder_save);
end

% loop over the frames
for frame_flag = 1:size(subImgs, 5)

    % make folder if it does not exist
    folder_write = fullfile(folder_save, ['Frame', num2str(frame_flag, '%04d')]);
    if ~isfolder(folder_write)
        mkdir(folder_write);
    end

    % loop over the subapertures
    for subaperture_flag = 1:size(subImgs, 4)
        file_write = fullfile(folder_write, ['Subaperture', num2str(subaperture_flag, '%02d'), '.mat']);

        % X is subaperture TORE volume, Y is proxy truth spot position
        X = squeeze(subImgs(:, :, :, subaperture_flag, frame_flag));
        Y = squeeze(positions(:, subaperture_flag, frame_flag))';

        % write the file with X, Y, dataset ID, subaperture ID, and frame ID
        save(file_write, 'X', 'Y', 'dataset_flag', 'subaperture_flag', 'frame_flag');
    end
end

fprintf(1, 'Done. '); toc;

%% plot tracks making sure they're okay after subsampling
% 
% fig1 = figure('Position', [100, 100, 1000, 400]);
% hold on;
% for ii = 1:numSubaps
% % for ii = 26
%     plot(tvec_ev, squeeze(positions(1,ii,:)), 'LineWidth', 2);
%     hold on;
%     plot(tvec_ev, squeeze(positions_baseline(1,ii,:)), 'LineWidth', 2);
%     xlim([tvec_ev(1), tvec_ev(end)]);
%     ylim([-25, 25]);
%     xlabel('Time');
%     ylabel('x-slope [event pixels]');
%     title(['X-Slope, Subaperture ', num2str(ii, '%02d')]);
%     drawnow;
% 
%     pause(0.75);
%     clf(fig1);
% end
% close(fig1);
% 
% fig2 = figure('Position', [100, 100, 1000, 400]);
% hold on;
% for ii = 1:numSubaps
% % for ii = 29
%     plot(tvec_ev, squeeze(positions(2,ii,:)), 'LineWidth', 2);
%     hold on;
%     plot(tvec_ev, squeeze(positions_baseline(2,ii,:)), 'LineWidth', 2);
%     xlim([tvec_ev(1), tvec_ev(end)]);
%     ylim([-25, 25]);
%     xlabel('Time');
%     ylabel('y-slope [event pixels]');
%     title(['Y-Slope, Subaperture ', num2str(ii, '%02d')]);
%     drawnow;
% 
%     pause(0.75);
%     clf(fig2);
% end
% close(fig2);
