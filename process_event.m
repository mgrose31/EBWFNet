clear; close all; clc;
addpath('..\..\Filtering');
addpath('C:\Users\ISSL\Documents\ISSL Sync\Mitchell\Scripts');

% Script that iterates over datasets (sequences) and builds TORE volumes at
% specified times (positive trigger events corresponding to FLIR frames).
% TORE volumes and meta data are saved for further processing.
% Author: Mitchell Grose, University of Dayton, 2023

%% script hyperparameters

evFrames_FPS = 100; % desired event frames FPS (want one-to-one match)
FLIR_FPS = 100; % FLIR frames captured at this FPS

subap_sz_ev = 54; % Prophesee VGA sensor subaperture sizes

pitch_EV = 15e-6; % Prophesee VGA sensor pixel pitch, 15 micrometers

%% data to use

date_flag = '2023-07-30';

% % loop over datasets (sequences)
for ii = 35
% for ii = [1:11, 13:35]

    % current dataset (sequence)
    dataset_flag = ii;

    fprintf(1, ['Dataset ', num2str(dataset_flag, '%02d'), '...\n']);

    % full filepath to dataset (sequence)
    fpath_ev = ['D:\Data\Dual_WFS\', date_flag, '_Eastwood\DualSync\Event\dataset', num2str(dataset_flag, '%02d')];

    % load the event data that has been saved as a .mat
    load(fullfile(fpath_ev, 'eventData_cd.mat'));

    % first frame to keep (defaults to 1)
    firstframe_keep = 1;

    %% process event data

    % add 1 to events due to MATLAB indexing
    events.x = events.x + 1;
    events.y = events.y + 1;

    % keep events within trigger events +/- 100 ms
    idx = events.ts >= events.trigger.ts(1) - 100000 & events.ts <= events.trigger.ts(end) + 100000;
    events.x = events.x(idx);
    events.y = events.y(idx);
    events.p = events.p(idx);
    events.ts = events.ts(idx);

    % remove hot pixels
    flag_RemoveHotPixels = true;
    if flag_RemoveHotPixels
        tic; fprintf(1, 'Removing hot pixels... ');

        % load hot pixel positions
        load('positions_hot_pixels.mat');

        % find events to remove
        idx_hot = zeros(size(events.ts));
        for jj=1:length(positions_hot_pixel)
            idx_hot_tmp = events.x == positions_hot_pixel(jj,1) & events.y == positions_hot_pixel(jj,2);
            idx_hot = idx_hot | idx_hot_tmp;
        end

        clear idx_hot_tmp;

        events.ts = events.ts(~idx_hot);
        events.p = events.p(~idx_hot);
        events.x = events.x(~idx_hot);
        events.y = events.y(~idx_hot);

        fprintf(1, 'Done! '); toc;
    end

    % Inceptive Event (IE) filter
    % false is for "no IE" data
    flag_IE_filter = true;
    if flag_IE_filter
        tic;

        % twindow_IE_filter = 2500; % 2.5 ms window
        twindow_IE_filter = 5000; % 5 ms window
        % twindow_IE_filter = 10000; % 10 ms window

        fprintf(1, ['Event filtering with ', num2str(twindow_IE_filter), ' microsecond window... ']);

        [IE, IEm] = IE_filter(events, twindow_IE_filter);

        %     IE = IE & (IEm > 1);

        % only keep IE events
        events.x = events.x(IE);
        events.y = events.y(IE);
        events.p = events.p(IE);
        events.ts = events.ts(IE);
        events.IEm = IEm(IE);

        fprintf(1,'Done! '); toc;
    else
        disp('Not using IE filter!');
        twindow_IE_filter = [];
    end

    %% make plots

    % accumulate all events
    evFrame_ref = accumarray([events.y, events.x], 1, [events.height, events.width], @sum, 0);

    % show event frame with accumulated events
    figure;
    imagesc(evFrame_ref);
    axis image;
    colormap(gray);
    colorbar;
    drawnow;

    % % plots for analysis of data processing
    % % idx_sub = ones(size(events.ts));
    % idx_sub = events.ts >= 1e6 & events.ts <= 2e6;
    % % idx_sub = events.ts >= events.trigger.ts(1) & events.ts <= events.trigger.ts(end);
    % idx_sub_pos = idx_sub & events.p == 1;
    % idx_sub_neg = idx_sub & events.p == -1;
    % 
    % figure;
    % scatter3(events.x(idx_sub_pos), events.y(idx_sub_pos), events.ts(idx_sub_pos), 1, 'blue');
    % hold on;
    % scatter3(events.x(idx_sub_neg), events.y(idx_sub_neg), events.ts(idx_sub_neg), 1, 'red');
    % % ylim([0, width]); zlim([0, height]);
    % % xlim([width/2 - 100, width/2 + 100]); ylim([height/2 - 100, height/2 + 100]);
    % xlim([75, 475]); ylim([0, 400]);
    % % view(45, 70);
    % title('Positive Event: Blue; Negative Event: Red');
    % xlabel('x');
    % ylabel('y');
    % zlabel('time (\mus)');
    % drawnow;
    % hold off;
    % 
    % % figure;
    % figure('Position', [100, 100, 1000, 400]);
    % scatter3(events.ts(idx_sub_pos), events.x(idx_sub_pos), events.y(idx_sub_pos), 1, 'blue');
    % hold on;
    % scatter3(events.ts(idx_sub_neg), events.x(idx_sub_neg), events.y(idx_sub_neg), 1, 'red');
    % ylim([0, events.width]); zlim([0, events.height]);
    % % ylim([width/2 - 100, width/2 + 100]); zlim([height/2 - 100, height/2 + 100]);
    % title('Positive Event: Blue; Negative Event: Red');
    % xlabel('time (\mus)');
    % ylabel('x');
    % zlabel('y');
    % drawnow;
    % hold off;

    %% define event data time vector

    % set event data time vector to positive trigger events
    tvec_ev = events.trigger.ts(events.trigger.p==1);

    % % add time for TORE volume generation
    toffset_ev = 3500; % 7500 theoretical, but 3500 works better
    tvec_ev = tvec_ev + toffset_ev;

    % height and width of event sensor
    height = events.height;
    width = events.width;

    %% build TORE volumes

    tic;
    evFrames_TORE = events2ToreFeature(events.x, events.y, events.ts, events.p, tvec_ev, 4, [height, width]);
    fprintf(1, 'Done. '); toc;

    %% save TORE volumes
    % uncomment code below to actually save

    tic; fprintf(1, 'Saving... ');
    save(['D:\tracks_ev_dataset', num2str(dataset_flag, '%02d'), '.mat'],...
        'evFrames_TORE', 'events', 'tvec_ev', 'toffset_ev', 'flag_IE_filter',...
        'twindow_IE_filter', 'date_flag', 'dataset_flag', 'evFrames_FPS', '-v7.3');
    fprintf(1, 'Done. '); toc;
    
    %% clear variables and close figures

    clear('evFrames_TORE', 'events', 'tvec_ev', 'toffset_ev',...
        'flag_IE_filter', 'twindow_IE_filter', 'dataset_flag');

    close all;

end