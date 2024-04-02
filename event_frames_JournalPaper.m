clear; close all; clc;

%% function paths
addpath('../Loading/');
addpath('../Filtering/');

%% set up filename of data to load

dataset_flag = 4;

% fpath = ['D:\Data\Dual_WFS\2023-07-30_Eastwood\DualSync\Event\dataset', num2str(dataset_flag, '%02d')];
% fname = 'EventData_cd.mat';
% 
% filename = fullfile(fpath, fname);
% 
% % event camera subaperture size
% subap_sz_ev = 50;
% 
% % FLIR camera subaperture size
% subap_sz_flir = 225;
% 
% %% load events
% 
% tic; fprintf(1, 'Loading events... ');
% load(filename);
% fprintf(1, 'Done! '); toc;
% 
% %% sensor dimensions
% % width = 640; % VGA sensor width
% % height = 480; % VGA sensor height
% 
% % width = 1280; % HD sensor width
% % height = 720; % HD sensor height
% 
% % use when going from .mat file saved by "save_events_as_mat.m"
% width = events.width;
% height = events.height;
% 
% %% start at time 0 and index events for Matlab
% 
% idx = events.ts >= events.trigger.ts(1) & events.ts <= events.trigger.ts(end);
% events.ts = events.ts(idx);
% events.x = events.x(idx);
% events.y = events.y(idx);
% events.p = events.p(idx);
% 
% % set starting time to 0 microseconds
% events.ts = events.ts - events.ts(1);
% 
% % add 1 to events due to Matlab indexing
% events.x = events.x + 1;
% events.y = events.y + 1;
% 
% %% downselect data
% 
% % % only use a portion of the sensor
% % idx = events.x >= 50 & events.x <= 450 ...
% %     & events.y >= 1 & events.y <= 400;
% % events.ts = events.ts(idx);
% % events.x = events.x(idx);
% % events.y = events.y(idx);
% % events.p = events.p(idx);
% 
% %% filter data (hot pixels and IE)
% 
% remove_hotPixels = 1;
% if remove_hotPixels
%     tic; fprintf(1, 'Removing hot pixels... ');
% 
%     load('../Tracking/WFS_Spots/positions_hot_pixels.mat');
% 
%     idx_hot = zeros(size(events.ts));
% 
%     for ii=1:length(positions_hot_pixel)
%         idx_hot_tmp = events.x == positions_hot_pixel(ii,1) & events.y == positions_hot_pixel(ii,2);
%         idx_hot = idx_hot | idx_hot_tmp;
%     end
% 
%     events.x = events.x(~idx_hot);
%     events.y = events.y(~idx_hot);
%     events.p = events.p(~idx_hot);
%     events.ts = events.ts(~idx_hot);
% 
%     fprintf(1, 'Done! '); toc;
% end
% 
% use_IE_Filter = 1;
% if use_IE_Filter
%     tic; fprintf(1,'Event filtering... ');
% 
%     twindow = 5000;
%     [IE, IEm] = IE_filter(events, twindow);
% 
%     fprintf(1,'Done! '); toc;
% 
% %     IE = IE & (IEm > 1);
% 
%     events.x = events.x(IE);
%     events.y = events.y(IE);
%     events.p = events.p(IE);
%     events.ts = events.ts(IE);
% end
% 
% %% build frames from events
% 
% % desired effective frames per second
% evFrames_FPS = 100;
% 
% % time step interval and window
% tint = (1/evFrames_FPS)*1e6;
% twin = (10e-3)*1e6; % (+/- 5 ms window)
% 
% tvec = events.ts(1):tint:events.ts(end);
% 
% evFrames = zeros(height, width, length(tvec));
% 
% for ii = 1:length(tvec)
%     twin_tmp = [tvec(ii) - twin, tvec(ii)];
%     idx = (events.ts >= twin_tmp(1)) & (events.ts < twin_tmp(2));
%     B = accumarray([events.y(idx), events.x(idx)], events.p(idx), [height, width], @sum, 0);
% 
%     evFrames(:,:,ii) = B;
% 
%     clear('B');
% end
% 
% figure;
% imagesc(evFrames(:,:,1000));
% axis image;
% colormap(gray);
% colorbar;
% clim([-1,1]);

%% event camera frame

% load(['../Tracking/WFS_Spots/spot_ref_centers_20230730/ref_pos_ev_dataset', num2str(dataset_flag, '%02d'), '.mat']);
% 
% % rectangles to draw for each subaperture position
% r_pos_ev = [xc_ref_ev-subap_sz_ev/2, yc_ref_ev-subap_sz_ev/2,...
%     ones(size(yc_ref_ev)).*subap_sz_ev, ones(size(yc_ref_ev)).*subap_sz_ev];
% 
% fig_event = figure;
% imagesc(B);
% colormap('gray');
% cbh = colorbar;
% cbh.Ticks = [-1,0,1];
% cbh.TickLabels = num2cell([-1,0,1]);
% clim([-1, 1]);
% axis image off xy;
% set(gca, 'FontWeight', 'bold');
% hold on;
% for ii=1:length(xc_ref_ev)
%     % rectangle('Position', r_pos_ev(ii,:), 'EdgeColor', 'g', 'LineWidth', 3);
%     rectangle('Position', r_pos_ev2(ii,:), 'EdgeColor', 'g', 'LineWidth', 2);
% end
% hold off;
% drawnow;
% 
% % exportgraphics(fig_event, 'Event_Frame_370_Dataset04.pdf');

%%

% Need to use tvec from formatted (?) file.
% Need to remove frames I throw out so I can put the target "x"

load(['C:\Users\ISSL\Documents\ISSL Sync\Mitchell\Scripts\MATLAB Prophesee\Tracking\WFS_Spots\Formatted_Data\formatted_dataset', num2str(dataset_flag, '%02d'), '.mat'], ...
    'events', 'tvec_ev', 'xc_ref_ev_rounded', 'yc_ref_ev_rounded', 'positions');

positions_Frame = positions;
clear('positions');

height = events.height;
width = events.width;

%% build frames from events

% desired effective frames per second
evFrames_FPS = 100;

% time step interval and window
tint = (1/evFrames_FPS)*1e6;
twin = (10e-3)*1e6; % (+/- 5 ms window)

evFrames = zeros(height, width, length(tvec_ev));

for ii = 1:length(tvec_ev)
    twin_tmp = [tvec_ev(ii) - twin, tvec_ev(ii)];
    idx = (events.ts >= twin_tmp(1)) & (events.ts < twin_tmp(2));
    B = accumarray([events.y(idx), events.x(idx)], events.p(idx), [height, width], @sum, 0);

    evFrames(:,:,ii) = B;

    clear('B');
end

figure;
imagesc(evFrames(:,:,1000));
axis image;
colormap(gray);
colorbar;
clim([-1,1]);

%%

[xx, yy] = meshgrid(-25:25-1, -25:25-1); % define sub-aperture size

NumFrames = size(evFrames, 3);

% subImgs = NaN([size(xx), size(evFrames_TORE, 3), length(xc_ref_ev_rounded), NumFrames]);
subImgs = NaN([size(xx), length(xc_ref_ev_rounded), NumFrames]);

% positions = NaN([2, length(xc_ref_ev), NumFrames]);

cnt = 0;
for ii = 1:NumFrames
    fprintf(1, ['Processing frame ', num2str(ii), ' of ', num2str(NumFrames), '\n']);
    for jj = 1:length(xc_ref_ev_rounded)

        img_tmp = double(evFrames(yc_ref_ev_rounded(jj)-25:yc_ref_ev_rounded(jj)+25-1,...
            xc_ref_ev_rounded(jj)-25:xc_ref_ev_rounded(jj)+25-1,ii));

        subImgs(:,:,jj,ii) = img_tmp;

        % positions(1,jj,ii) = xc_imgs_flirT_new(jj,ii);
        % positions(2,jj,ii) = yc_imgs_flirT_new(jj,ii);

    end
end

%%

idx_plot_frame = 300;
idx_plot_subap = 8;

figure;
imagesc(xx(1,:), yy(:,1), subImgs(:, :, idx_plot_subap, idx_plot_frame));
axis image xy;
colormap(gray); cb = colorbar; clim([-1, 1]);
cb.Ticks = -1:1;
hold on;
plot(positions_Frame(1, idx_plot_subap, idx_plot_frame), positions_Frame(2, idx_plot_subap, idx_plot_frame), 'rx', 'LineWidth', 3);
set(gca, 'FontWeight', 'Bold');

%%

fig1 = figure;
set(gca, 'FontWeight', 'Bold');
hold on;
for ii = 1:50
    idx_plot_frame = ii;
    imagesc(xx(1,:), yy(:,1), subImgs(:, :, idx_plot_subap, idx_plot_frame));
    axis image xy;
    colormap(gray); cb = colorbar; clim([-1, 1]);
    cb.Ticks = -1:1;
    hold on;
    plot(positions_Frame(1, idx_plot_subap, idx_plot_frame), positions_Frame(2, idx_plot_subap, idx_plot_frame), 'rx', 'LineWidth', 3);
    hold off;

    drawnow;
    pause(0.50);

    clf(fig1);
end
close(fig1);

%%

idx_plot_subap = ones(1, 16) .* 9;
% idx_plot_frame = 300:300 + 15;
idx_plot_frame = 28:40:size(subImgs, 4);
idx_plot_frame = idx_plot_frame(1:16);

% idx_plot_subap = 8:8+15;
% idx_plot_frame = ones(1, 16) .* 200;


fig2 = figure('Position', [300, 300, 800, 600]);
% tiledlayout(4, 4);
tiledlayout(3, 4, 'TileSpacing', 'compact');
% for ii = 1:16
for ii = 1:12
    nexttile;

    imagesc(xx(1,:), yy(:,1), subImgs(:, :, idx_plot_subap(ii), idx_plot_frame(ii)));
    axis image xy off;
    colormap(gray);
    % cb = colorbar; clim([-1, 1]); cb.Ticks = -1:1;
    clim([-1, 1]);
    hold on;
    plot(positions_Frame(1, idx_plot_subap(ii), idx_plot_frame(ii)), ...
        positions_Frame(2, idx_plot_subap(ii), idx_plot_frame(ii)), 'rx', 'LineWidth', 3);
    set(gca, 'FontWeight', 'Bold');
    hold off;
end

% exportgraphics(fig2, 'Example_Subaperture_Event_Frames.pdf');

%%

fig3 = figure;
% for ii = 1:12
for ii = 12
    imagesc(xx(1,:), yy(:,1), subImgs(:, :, idx_plot_subap(ii), idx_plot_frame(ii)));
    axis image xy off;
    colormap(gray);
    % cb = colorbar; clim([-1, 1]); cb.Ticks = -1:1;
    clim([-1, 1]);
    hold on;
    plot(positions_Frame(1, idx_plot_subap(ii), idx_plot_frame(ii)), ...
        positions_Frame(2, idx_plot_subap(ii), idx_plot_frame(ii)), 'rx', 'LineWidth', 3, 'MarkerSize', 15);
    set(gca, 'FontWeight', 'Bold');
    hold off;

    % exportgraphics(fig3, ['Figure_Example_Subaperture_Event_Frame', num2str(ii, '%02d'), '.pdf']);

end