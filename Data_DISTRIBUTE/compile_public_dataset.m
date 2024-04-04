clear; close all; clc;
addpath('../WFS_Spots/');
addpath('../../Filtering');

%%

dataset_flag_arr = [1:8, 10, 11, 13:17, 19, 22, 24:35];

num_datasets = length(dataset_flag_arr);

wb = waitbar(0.);

cnt = 0;
for dataset_flag = dataset_flag_arr(1:end)
    cnt = cnt + 1;
    waitbar(cnt/num_datasets, wb, ['Processing dataset ', num2str(cnt), ' of ', num2str(num_datasets)]);
    
    process_data(dataset_flag);
end
close(wb);


function process_data(dataset_flag)

fpath_events = ['D:\Data\Dual_WFS\2023-07-30_Eastwood\DualSync\Event\dataset', num2str(dataset_flag, '%02d'), '\EventData_cd.mat'];

fpath_formatted = ['C:\Users\ISSL\Documents\ISSL Sync\Mitchell\Scripts\MATLAB Prophesee\Tracking\WFS_Spots\Formatted_Data\', ...
    'formatted_dataset', num2str(dataset_flag, '%02d'), '.mat'];

% get the original event stream
load(fpath_events, 'events');

% get the information from data formatting
load(fpath_formatted, 'xc_ref_ev_rounded', 'yc_ref_ev_rounded', 'positions', 'positions_baseline', 'tvec_ev', 'frameID');

xPos_Frame_tmp = squeeze(positions(1,:,:));
yPos_Frame_tmp = squeeze(positions(2,:,:));

xPos_AR1_tmp = squeeze(positions_baseline(1,:,:));
yPos_AR1_tmp = squeeze(positions_baseline(2,:,:));

clear('positions', 'positions_baseline');

[numSubaps, numFrames] = size(xPos_Frame_tmp);

%% process event data

% % add 1 to events due to MATLAB indexing
events.x = events.x + 1;
events.y = events.y + 1;

% keep events within trigger events +/- 1 second
idx_frame = events.ts >= events.trigger.ts(1) - 1e6 ...
    & events.ts <= events.trigger.ts(end) + 1e6;
events.x = events.x(idx_frame);
events.y = events.y(idx_frame);
events.p = events.p(idx_frame);
events.ts = events.ts(idx_frame);

% add time to be better aligned with frame data
events.ts = events.ts + 3500;
events.trigger.ts = events.trigger.ts + 3500;

frames_keep = frameID;

events_original = events;
clear('events');

%%

% loop over subaperture
for idx_subap = 1:numSubaps
    fprintf(1, ['Processing subaperture ', num2str(idx_subap), ' of ', num2str(numSubaps), '.\n']);
    % get events within subaperture
    idx_events_subap = events_original.y >= yc_ref_ev_rounded(idx_subap) - 25 ...
        & events_original.y <= yc_ref_ev_rounded(idx_subap) + 24 ...
        & events_original.x >= xc_ref_ev_rounded(idx_subap) - 25 ...
        & events_original.x <= xc_ref_ev_rounded(idx_subap) + 24;
    if sum(idx_events_subap) == 0
        error("No data in subaperture!");
    end

    % event stream for current subaperture
    events = events_original;

    subap_center_x = xc_ref_ev_rounded(idx_subap);
    subap_center_y = yc_ref_ev_rounded(idx_subap);

    % only keep events in the current subaperture
    % center the events in the subaperture
    events.x = events.x(idx_events_subap) - subap_center_x;
    events.y = events.y(idx_events_subap) - subap_center_y;
    events.p = events.p(idx_events_subap);
    events.ts = events.ts(idx_events_subap);

    xPos_Frame = xPos_Frame_tmp(idx_subap,:);
    yPos_Frame = yPos_Frame_tmp(idx_subap,:);
    xPos_AR1 = xPos_AR1_tmp(idx_subap,:);
    yPos_AR1 = yPos_AR1_tmp(idx_subap,:);

    % figure;
    % plot(tvec_ev, xPos_Frame, 'LineWidth', 2);
    % hold on;
    % plot(tvec_ev, xPos_AR1, 'LineWidth', 2);

    subaperture_flag = idx_subap;

    folder_save = ['./Data_DISTRIBUTE/Dataset', num2str(dataset_flag, '%02d')];
    if ~isfolder(folder_save)
        mkdir(folder_save);
    end

    save([folder_save, '/Formatted_Dataset', num2str(dataset_flag, '%02d'), ...
        '_Subaperture', num2str(idx_subap, '%02d'), '_DISTRIBUTE.mat'], ...
    'events', 'dataset_flag', 'frames_keep', 'numFrames', ...
    'tvec_ev', 'subap_center_x', 'subap_center_y', 'subaperture_flag', ...
    'xPos_AR1', 'yPos_AR1', 'xPos_Frame', 'yPos_Frame');


end



end
