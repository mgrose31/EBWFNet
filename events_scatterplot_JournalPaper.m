clear; close all; clc;

%% function paths
addpath('../Loading/');
addpath('../Filtering/');

%% set up filename of data to load

fpath = 'D:\Data\Dual_WFS\2023-07-30_Eastwood\DualSync\Event\dataset04';
fname = 'EventData_cd.mat';

filename = fullfile(fpath, fname);

%% load events

tic; fprintf(1, 'Loading events... ');
% events = load_cd_events(filename, false, false, Inf);
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

use_IE_Filter = 1;
if use_IE_Filter
    tic; fprintf(1,'Event filtering... ');

    % twindow= 10000;
    twindow = 5000;
    %IE = findInceptiveEvents(events, twindow);
    [IE, IEm] = IE_filter(events, twindow);

    fprintf(1,'Done! '); toc;

%     IE = IE & (IEm > 1);

    events.x = events.x(IE);
    events.y = events.y(IE);
    events.p = events.p(IE);
    events.ts = events.ts(IE);
end

%% make scatter plots

% idx_sub = ones(size(events.ts));
% idx_sub = events.ts >= 5e6 & events.ts <= 6e6;
% idx_sub = events.ts >= 5e6 & events.ts <= 5.5e6;
idx_sub = events.ts >= 5e6 & events.ts <= 5.1e6;
% idx_sub = events.ts >= events.trigger.ts(1) & events.ts <= events.trigger.ts(end);
idx_sub_pos = idx_sub & events.p == 1;
idx_sub_neg = idx_sub & events.p == -1;

fig1 = figure;
scatter3(events.x(idx_sub_pos), events.y(idx_sub_pos), events.ts(idx_sub_pos) ./ 1e6, 1, 'blue');
hold on;
scatter3(events.x(idx_sub_neg), events.y(idx_sub_neg), events.ts(idx_sub_neg) ./ 1e6, 1, 'red');
% xlim([0, width]); ylim([0, height]);
xlim([235, 325]);
ylim([50, 175]);
view(-37, 30);
xlabel('x (pixels)');
ylabel('y (pixels)');
zlabel('time (s)');
xticks(240:20:320);
yticks(60:20:160);
zticks(5:0.02:5.1);
grid on;
set(gca, 'FontWeight', 'bold', 'FontSize', 12);
hold off;

% exportgraphics(fig1, 'EventData1.pdf');
% exportgraphics(fig1, 'EventData_Raw.pdf');
% exportgraphics(fig1, 'EventData_HotPixelRemoved.pdf');
% exportgraphics(fig1, 'EventData_HotPixelRemoved_IE_05ms.pdf');


[X_test, Y_test] = meshgrid(1:width, 1:height);
% Z_test = ones(size(X_test)).*5;
% Z_test2 = ones(size(X_test)).*5 + 0.001;


fig2 = figure;
scatter3(events.x(idx_sub_pos), events.y(idx_sub_pos), events.ts(idx_sub_pos) ./ 1e6, 1, 'blue');
hold on;
scatter3(events.x(idx_sub_neg), events.y(idx_sub_neg), events.ts(idx_sub_neg) ./ 1e6, 1, 'red');
for ii = 1:10
    Z_test = ones(size(X_test)) .* 5 + (ii-1) * 0.01;
    surf(X_test, Y_test, Z_test, 'LineWidth', 1.5);
end
view(0, 0);
xlim([235, 325]);
ylim([50, 110]);
xlabel('x (pixels)');
ylabel('y (pixels)');
zlabel('time (s)');
xticks(240:10:320);
yticks(60:20:160);
zticks(5:0.01:5.1);
set(gca, 'FontWeight', 'bold', 'FontSize', 12);
hold off;

% exportgraphics(fig2, 'EventData2.pdf');
% exportgraphics(fig2, 'EventData_HotPixelRemoved_IE_05ms_Frames_100Hz.pdf');


figure;
scatter3(events.x(idx_sub_pos), events.y(idx_sub_pos), events.ts(idx_sub_pos) ./ 1e6, 1, 'blue');
hold on;
scatter3(events.x(idx_sub_neg), events.y(idx_sub_neg), events.ts(idx_sub_neg) ./ 1e6, 1, 'red');
xlim([0, width]); ylim([0, height]);
view(0, 90);
xlabel('x (pixels)');
ylabel('y (pixels)');
zlabel('time (s)');
set(gca, 'FontWeight', 'bold');
hold off;

% % figure;
% figure('Position', [100, 100, 1200, 400]);
% scatter3(events.ts(idx_sub_pos), events.x(idx_sub_pos), events.y(idx_sub_pos), 1, 'blue');
% hold on;
% scatter3(events.ts(idx_sub_neg), events.x(idx_sub_neg), events.y(idx_sub_neg), 1, 'red');
% ylim([0, width]); zlim([0, height]);
% % ylim([width/2 - 100, width/2 + 100]); zlim([height/2 - 100, height/2 + 100]);
% % title('Positive Event: Blue; Negative Event: Red');
% xlabel('time (\mus)', 'FontWeight', 'bold');
% ylabel('x', 'FontWeight', 'bold');
% zlabel('y', 'FontWeight', 'bold');
% set(gca, 'FontWeight', 'bold');
% hold off;

%% look at the trigger events

% generate planes at the trigger frames
[X_test, Y_test] = meshgrid(1:500, 1:500);
Z_test = ones(size(X_test)).*events.trigger.ts(1);
Z_test2 = ones(size(X_test)).*events.trigger.ts(1) + 1000;
Z_test3 = ones(size(X_test)).*events.trigger.ts(3);
Z_test4 = ones(size(X_test)).*events.trigger.ts(3) + 1000;
Z_test5 = ones(size(X_test)).*events.trigger.ts(5);
Z_test6 = ones(size(X_test)).*events.trigger.ts(5) + 1000;
Z_test7 = ones(size(X_test)).*events.trigger.ts(7);
Z_test8 = ones(size(X_test)).*events.trigger.ts(7) + 1000;
Z_test9 = ones(size(X_test)).*events.trigger.ts(9);
Z_test10 = ones(size(X_test)).*events.trigger.ts(9) + 1000;
Z_test11 = ones(size(X_test)).*events.trigger.ts(11);
Z_test12 = ones(size(X_test)).*events.trigger.ts(11) + 1000;
% Z_test2 = ones(size(X_test)).*events.trigger.ts(3);
% Z_test3 = ones(size(X_test)).*events.trigger.ts(5);
% Z_test4 = ones(size(X_test)).*events.trigger.ts(7);
% Z_test5 = ones(size(X_test)).*events.trigger.ts(9);
% Z_test6 = ones(size(X_test)).*events.trigger.ts(11);

% idx_sub = ones(size(events.ts));
idx_sub = events.ts >= events.trigger.ts(1) ...
    & events.ts < events.trigger.ts(11);
idx_sub_pos = idx_sub & events.p == 1;
idx_sub_neg = idx_sub & events.p == -1;

figure;
scatter3(events.x(idx_sub_pos), events.y(idx_sub_pos), events.ts(idx_sub_pos), 1, 'blue');
hold on;
scatter3(events.x(idx_sub_neg), events.y(idx_sub_neg), events.ts(idx_sub_neg), 1, 'red');
surf(X_test, Y_test, Z_test, 'LineWidth', 1.5);
surf(X_test, Y_test, Z_test2, 'LineWidth', 1.5);
surf(X_test, Y_test, Z_test3, 'LineWidth', 1.5);
surf(X_test, Y_test, Z_test4, 'LineWidth', 1.5);
surf(X_test, Y_test, Z_test5, 'LineWidth', 1.5);
surf(X_test, Y_test, Z_test6, 'LineWidth', 1.5);
surf(X_test, Y_test, Z_test7, 'LineWidth', 1.5);
surf(X_test, Y_test, Z_test8, 'LineWidth', 1.5);
surf(X_test, Y_test, Z_test9, 'LineWidth', 1.5);
surf(X_test, Y_test, Z_test10, 'LineWidth', 1.5);
surf(X_test, Y_test, Z_test11, 'LineWidth', 1.5);
surf(X_test, Y_test, Z_test12, 'LineWidth', 1.5);
% ylim([0, width]); zlim([0, height]);
% xlim([width/2 - 100, width/2 + 100]); ylim([height/2 - 100, height/2 + 100]);
xlim([0, width]); ylim([0, height]);
view(90, 0);
% title('Positive Event: Blue; Negative Event: Red');
xlabel('x');
ylabel('y');
zlabel('time (\mus)');
hold off;


%% make a movie of scatter3

% close all;
% 
% tint = 0.01e6;
% twindow = 0.2e6;
% t = 0;
% % tend = max(events.ts) + twindow;
% tend = 1e6;
% 
% v = VideoWriter('test', 'MPEG-4');
% v.Quality = 100;
% 
% v.FrameRate = 10;
% 
% open(v);
% 
% fig = figure(50);
% set(fig, 'Position', [100, 100, width, height]);
% hold on;
% 
% while t < tend
% 
%     idx_sub = (events.ts >= t - twindow / 2) & (events.ts < t + twindow / 2);
% %     idx_sub = ones(size(events.ts));
%     idx_sub_pos = idx_sub & events.p == 1;
%     idx_sub_neg = idx_sub & events.p == -1;
% 
% %     figure('Position',[100,100,width,height]);
%     % scatter3(events.ts(idx_sub), events.x(idx_sub), events.y(idx_sub), 1);
%     scatter3(events.x(idx_sub_pos), events.y(idx_sub_pos), events.ts(idx_sub_pos), 1, 'blue');
%     hold on;
%     scatter3(events.x(idx_sub_neg), events.y(idx_sub_neg), events.ts(idx_sub_neg), 1, 'red');
%     % ylim([0, width]); zlim([0, height]);
% %     xlim([width/2 - 100, width/2 + 100]); ylim([height/2 - 100, height/2 + 100]);
%     xlim([0, 500]);
%     ylim([0, 500]);
%     zlim([t - twindow / 2, t + twindow / 2]);
%     xlabel('x');
%     ylabel('y');
%     zlabel('time (\mus)');
% 
%     view(45, 70);
% 
%     writeVideo(v, getframe(gcf));
% 
% %     close all;
% 
%     t = t + tint;
% 
%     clf(fig);
% 
% end
% close(v);