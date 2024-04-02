clear; close all; clc;

%% function paths
addpath('../Loading/');
addpath('../Filtering/');

%% set up filename of data to load

dataset_flag = 4;

fpath = ['D:\Data\Dual_WFS\2023-07-30_Eastwood\DualSync\Event\dataset', num2str(dataset_flag, '%02d')];
fname = 'EventData_cd.mat';

filename = fullfile(fpath, fname);

% event camera subaperture size
subap_sz_ev = 50;

% FLIR camera subaperture size
subap_sz_flir = 225;

%% load events

tic; fprintf(1, 'Loading events... ');
load(filename);
fprintf(1, 'Done! '); toc;

%% sensor dimensions
% width = 640; % VGA sensor width
% height = 480; % VGA sensor height

% width = 1280; % HD sensor width
% height = 720; % HD sensor height

% use when going from .mat file saved by "save_events_as_mat.m"
width = events.width;
height = events.height;

%% start at time 0 and index events for Matlab

idx = events.ts >= events.trigger.ts(1) & events.ts <= events.trigger.ts(end);
events.ts = events.ts(idx);
events.x = events.x(idx);
events.y = events.y(idx);
events.p = events.p(idx);

% set starting time to 0 microseconds
events.ts = events.ts - events.ts(1);

% add 1 to events due to Matlab indexing
events.x = events.x + 1;
events.y = events.y + 1;

%% downselect data

% % only use a portion of the sensor
% idx = events.x >= 50 & events.x <= 450 ...
%     & events.y >= 1 & events.y <= 400;
% events.ts = events.ts(idx);
% events.x = events.x(idx);
% events.y = events.y(idx);
% events.p = events.p(idx);

%% filter data (hot pixels and IE)

remove_hotPixels = 1;
if remove_hotPixels
    tic; fprintf(1, 'Removing hot pixels... ');

    load('../Tracking/WFS_Spots/positions_hot_pixels.mat');

    idx_hot = zeros(size(events.ts));

    for ii=1:length(positions_hot_pixel)
        idx_hot_tmp = events.x == positions_hot_pixel(ii,1) & events.y == positions_hot_pixel(ii,2);
        idx_hot = idx_hot | idx_hot_tmp;
    end

    events.x = events.x(~idx_hot);
    events.y = events.y(~idx_hot);
    events.p = events.p(~idx_hot);
    events.ts = events.ts(~idx_hot);

    fprintf(1, 'Done! '); toc;
end

use_IE_Filter = 0;
if use_IE_Filter
    tic; fprintf(1,'Event filtering... ');

    twindow = 5000;
    [IE, IEm] = IE_filter(events, twindow);

    fprintf(1,'Done! '); toc;

%     IE = IE & (IEm > 1);

    events.x = events.x(IE);
    events.y = events.y(IE);
    events.p = events.p(IE);
    events.ts = events.ts(IE);
end

%% load FLIR data

thresh_global = 0.3;

fpath_flir = ['D:\Data\Dual_WFS\2023-07-30_Eastwood\DualSync\FLIR\dataset', num2str(dataset_flag, '%02d')];

dd_FLIR = dir(fullfile(fpath_flir, '*.tiff'));

imgs = NaN(1536, 2048, 500);

tic; fprintf(1, 'Reading FLIR images... ');
for ii=1:size(imgs, 3)
    imgs(:,:,ii) = fliplr(double(imread(fullfile(dd_FLIR(ii).folder, dd_FLIR(ii).name))));
end
fprintf(1, 'Done. '); toc;


% apply global threshold to images
tic; fprintf(1, 'Filtering images with global threshold... ');
imgs_thresh = imgs;
clear('imgs');

for ii=1:size(imgs_thresh, 3)
    imgs_thresh(:,:,ii) = imgs_thresh(:,:,ii) - max(imgs_thresh(:,:,ii), [], 'all').*thresh_global;
end
imgs_thresh(imgs_thresh<0) = 0;
fprintf(1, 'Done. '); toc;


%% get data to visualize

% frame of data we want
frame_plt = 370;

% find indices corresponding to positive trigger events
idx_trigger_pos = find(events.trigger.p==1);

% get the timestamp of the "frame_plt" positive trigger event
ts_trigger = events.trigger.ts(idx_trigger_pos(frame_plt));

% get data 10ms before and 1ms after the trigger; add 3500s for offset
idx_sub = events.ts >= (ts_trigger + 3500) - 1e4 & events.ts <= (ts_trigger + 3500) + 1e3;

% show the FLIR frame
fig_frame = figure;
imagesc(imgs_thresh(:,:,frame_plt));
axis image xy off;
colormap(gray);

% build the event frame
B = accumarray([events.y(idx_sub), events.x(idx_sub)], events.p(idx_sub), [height, width], @sum, 0);

fig_event = figure;
imagesc(B);
colormap('gray');
% colorbar;
clim([-1, 1]);
axis image off xy

%% subaperture indices corresponding to each row and column

idx_row1 = [9, 16, 20, 23, 27];
idx_row2 = [4, 10, 17, 21, 24, 28, 33];
idx_row3 = [1, 5, 11, 29, 34, 38];
idx_row4 = [2, 6, 12, 35, 39];
idx_row5 = [3, 7, 13, 30, 36, 40];
idx_row6 = [8, 14, 18, 25, 31, 37];
idx_row7 = [15, 19, 22, 26, 32];

idx_col1 = 1:3;
idx_col2 = 4:8;
idx_col3 = 9:15;
idx_col4 = 16:19;
idx_col5 = 20:22;
idx_col6 = 23:26;
idx_col7 = 27:32;
idx_col8 = 33:37;
idx_col9 = 38:40;

%% load reference positions

load(['../Tracking/WFS_Spots/spot_ref_centers_20230730/ref_pos_flir_dataset', num2str(dataset_flag, '%02d'), '.mat']);

% define rectangles for subaperture positions
r_pos_flir = [xc_ref_flir-subap_sz_flir/2, yc_ref_flir-subap_sz_flir/2,...
    ones(size(yc_ref_flir)).*subap_sz_flir, ones(size(yc_ref_flir)).*subap_sz_flir];

% number of pixels to pad
pad_amt = 60;

% pad the images just for visualization of the subaperture boxes
flirFrame_ref_padded = padarray(imgs_thresh(:,:,frame_plt), [pad_amt, pad_amt], 0, 'both');
r_pos_flir_padded = r_pos_flir;
r_pos_flir_padded(:,1:2) = r_pos_flir_padded(:,1:2) + pad_amt;

r_pos_flir_padded2 = r_pos_flir_padded;
r_pos_flir_padded2(idx_row1, 2) = mean(r_pos_flir_padded(idx_row1, 2));
r_pos_flir_padded2(idx_row2, 2) = mean(r_pos_flir_padded(idx_row2, 2));
r_pos_flir_padded2(idx_row3, 2) = mean(r_pos_flir_padded(idx_row3, 2));
r_pos_flir_padded2(idx_row4, 2) = mean(r_pos_flir_padded(idx_row4, 2));
r_pos_flir_padded2(idx_row5, 2) = mean(r_pos_flir_padded(idx_row5, 2));
r_pos_flir_padded2(idx_row6, 2) = mean(r_pos_flir_padded(idx_row6, 2));
r_pos_flir_padded2(idx_row7, 2) = mean(r_pos_flir_padded(idx_row7, 2));
r_pos_flir_padded2(idx_col1, 1) = mean(r_pos_flir_padded(idx_col1, 1));
r_pos_flir_padded2(idx_col2, 1) = mean(r_pos_flir_padded(idx_col2, 1));
r_pos_flir_padded2(idx_col3, 1) = mean(r_pos_flir_padded(idx_col3, 1));
r_pos_flir_padded2(idx_col4, 1) = mean(r_pos_flir_padded(idx_col4, 1));
r_pos_flir_padded2(idx_col5, 1) = mean(r_pos_flir_padded(idx_col5, 1));
r_pos_flir_padded2(idx_col6, 1) = mean(r_pos_flir_padded(idx_col6, 1));
r_pos_flir_padded2(idx_col7, 1) = mean(r_pos_flir_padded(idx_col7, 1));
r_pos_flir_padded2(idx_col8, 1) = mean(r_pos_flir_padded(idx_col8, 1));
r_pos_flir_padded2(idx_col9, 1) = mean(r_pos_flir_padded(idx_col9, 1));

% r_pos_flir_padded2(idx_row1, 2) = r_pos_flir_padded2(idx_row1, 2) + 35;
% r_pos_flir_padded2(idx_row2, 2) = r_pos_flir_padded2(idx_row2, 2) + 15;
% r_pos_flir_padded2(idx_row7, 2) = r_pos_flir_padded2(idx_row7, 2) - 5;
% 
% r_pos_flir_padded2(idx_col1, 1) = r_pos_flir_padded2(idx_col1, 1) + 15;
% r_pos_flir_padded2(idx_col2, 1) = r_pos_flir_padded2(idx_col2, 1) + 15;

% plot reference image with subaperture positions
fig1 = figure;
% imagesc(imgs_thresh(:,:,frame_plt));
imagesc(flirFrame_ref_padded);
hold on;
% plot(xc_ref_flir, yc_ref_flir, 'r+', 'LineWidth', 2, 'MarkerSize', 6);
for ii=1:length(xc_ref_flir)
    % rectangle('Position', r_pos_flir(ii,:), 'EdgeColor', 'g', 'LineWidth', 2);
    % rectangle('Position', r_pos_flir_padded(ii,:), 'EdgeColor', 'g', 'LineWidth', 3);
    rectangle('Position', r_pos_flir_padded2(ii,:), 'EdgeColor', 'g', 'LineWidth', 2);
end
axis image xy off;
colorbar;
colormap('gray');
% title('Compute Subaperture Reference Positions');
set(gca, 'FontWeight', 'Bold');
hold off;
drawnow;

% exportgraphics(fig1, 'FLIR_Reference_Frame.pdf');
% exportgraphics(fig1, 'FLIR_Frame_370_Dataset04.pdf');


%% event camera frame

load(['../Tracking/WFS_Spots/spot_ref_centers_20230730/ref_pos_ev_dataset', num2str(dataset_flag, '%02d'), '.mat']);

% rectangles to draw for each subaperture position
r_pos_ev = [xc_ref_ev-subap_sz_ev/2, yc_ref_ev-subap_sz_ev/2,...
    ones(size(yc_ref_ev)).*subap_sz_ev, ones(size(yc_ref_ev)).*subap_sz_ev];

r_pos_ev2 = r_pos_ev;
r_pos_ev2(idx_row1, 2) = mean(r_pos_ev(idx_row1, 2));
r_pos_ev2(idx_row2, 2) = mean(r_pos_ev(idx_row2, 2));
r_pos_ev2(idx_row3, 2) = mean(r_pos_ev(idx_row3, 2));
r_pos_ev2(idx_row4, 2) = mean(r_pos_ev(idx_row4, 2));
r_pos_ev2(idx_row5, 2) = mean(r_pos_ev(idx_row5, 2));
r_pos_ev2(idx_row6, 2) = mean(r_pos_ev(idx_row6, 2));
r_pos_ev2(idx_row7, 2) = mean(r_pos_ev(idx_row7, 2));
r_pos_ev2(idx_col1, 1) = mean(r_pos_ev(idx_col1, 1));
r_pos_ev2(idx_col2, 1) = mean(r_pos_ev(idx_col2, 1));
r_pos_ev2(idx_col3, 1) = mean(r_pos_ev(idx_col3, 1));
r_pos_ev2(idx_col4, 1) = mean(r_pos_ev(idx_col4, 1));
r_pos_ev2(idx_col5, 1) = mean(r_pos_ev(idx_col5, 1));
r_pos_ev2(idx_col6, 1) = mean(r_pos_ev(idx_col6, 1));
r_pos_ev2(idx_col7, 1) = mean(r_pos_ev(idx_col7, 1));
r_pos_ev2(idx_col8, 1) = mean(r_pos_ev(idx_col8, 1));
r_pos_ev2(idx_col9, 1) = mean(r_pos_ev(idx_col9, 1));

fig_event = figure;
imagesc(B);
colormap('gray');
cbh = colorbar;
cbh.Ticks = [-1,0,1];
cbh.TickLabels = num2cell([-1,0,1]);
clim([-1, 1]);
axis image off xy;
set(gca, 'FontWeight', 'bold');
hold on;
for ii=1:length(xc_ref_ev)
    % rectangle('Position', r_pos_ev(ii,:), 'EdgeColor', 'g', 'LineWidth', 3);
    rectangle('Position', r_pos_ev2(ii,:), 'EdgeColor', 'g', 'LineWidth', 2);
end
hold off;
drawnow;

% exportgraphics(fig_event, 'Event_Frame_370_Dataset04.pdf');