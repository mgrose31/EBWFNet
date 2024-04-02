clear; close all; clc;
addpath('..\..\Filtering');

% Code to load event data (that has been written to a .mat file),
% accumulate events, then identify hot pixels and save the hot pixels as a
% .mat file for later filtering.
% Author: Mitchell Grose, University of Dayton, 2023

%% load all the data

date_flag = '2023-07-30';  % load data from this date

cnt = 0;
for ii = 1:35
    dataset_flag = ii;
    cnt = cnt + 1;

    fprintf(1, ['Dataset ', num2str(dataset_flag, '%02d'), '...\n']);

    fpath_ev = ['D:\Data\Dual_WFS\', date_flag, '_Eastwood\DualSync\Event\dataset', num2str(dataset_flag, '%02d')];

    % load event data
    dH5 = load(fullfile(fpath_ev, 'eventData_cd.mat'));
    
    if cnt == 1
        events = dH5.events;
    end

    % append event data
    events.ts = [events.ts; dH5.events.ts];
    events.x = [events.x; dH5.events.x];
    events.y = [events.y; dH5.events.y];
    events.p = [events.p; dH5.events.p];

end

%% accumulate events to visualize hot pixels

% for MATLAB indexing
events.x = events.x + 1;
events.y = events.y + 1;

evFrame_ref = accumarray([events.y, events.x], 1, [events.height, events.width], @sum, 0);

figure;
imagesc(evFrame_ref);
axis image;
colormap(gray);
colorbar;

%% 

% load prior hot pixels
load('positions_hot_pixels.mat');
% positions_hot_pixel = [positions_hot_pixel; vertcat(cursor_info.Position)];

% get CDF of accumulated events
[f1, x1] = ecdf(evFrame_ref(:));

% plot CDF of accumulated events
figure;
plot(x1, f1, 'LineWidth', 2);
xlabel('Number of Events');
ylabel('P[X \leq x]')

% create binary mask to visualize hot pixels
mask_off = ones(size(evFrame_ref));
for ii = 1:size(positions_hot_pixel, 1)
    mask_off(positions_hot_pixel(ii,2), positions_hot_pixel(ii,1)) = 0;
end

% set hot pixels to 0 in frame of accumulated events
evFrame_ref2 = evFrame_ref;
for ii = 1:size(positions_hot_pixel, 1)
    evFrame_ref2(positions_hot_pixel(ii,2), positions_hot_pixel(ii,1)) = 0;
end

% visualize frame of accumulated events with hot pixels set to 0
figure;
imagesc(evFrame_ref2);
axis image;
colormap(gray);
colorbar;

% get new CDF
[f1, x1] = ecdf(evFrame_ref2(:));

% plot new CDF
figure;
plot(x1, f1, 'LineWidth', 2);
xlabel('Number of Events');
ylabel('P[X \leq x]')

% find pixels with more than n events
[r, c] = find(evFrame_ref2 > 600000);
% [r, c] = find(evFrame_ref2 > 1000);

%% save hot pixels

% saving format is [x,y] which corresponds to [c,r]
positions_hot_pixel = [positions_hot_pixel; [c, r]]; % append hot pixels

% save hot pixels
% save('positions_hot_pixels.mat', 'positions_hot_pixel');
