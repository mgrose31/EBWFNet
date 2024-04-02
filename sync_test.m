clear; close all; clc;

%% function paths
% addpath('../Loading/');
addpath('../Filtering/');

%% set data paths

dd = dir('D:\Data\Sync_Tests\2023-07-06_TemporalDelay\FLIR\dataset01\*.tiff');

fpath_ev = 'D:\Data\Sync_Tests\2023-07-06_TemporalDelay\Event\dataset01';
fname_ev = 'EventData_cd.mat';

%% Set FLIR FPS

FLIR_FPS = 100;

% % uncomment if *not* using Read Out for Trigger Overlap Mode
% FLIR_FPS = 1 / (0.001 + (1 / FLIR_FPS));

tvec_FLIR = (0:length(dd)-1)./FLIR_FPS;

% framerate to generate event frames
evFrames_FPS = 500;


%% read the data

% imgs = NaN(1536, 2048, length(dd));
imgs = NaN(480, 640, length(dd));


tic; fprintf(1, 'Reading images...\n');
for ii=1:length(dd)
    fprintf(1, ['Reading image ', num2str(ii), ' of ', num2str(length(dd)), '\n']);

    % fliplr because 45 degree beamsplitter flips the image left-right
    imgs_tmp = fliplr(double(imread(fullfile(dd(ii).folder, dd(ii).name)))); 
    imgs(:,:,ii) = imresize(imgs_tmp, [480, 640]);
end
fprintf(1, 'Done. '); toc;

clear('imgs_tmp');

tic; fprintf(1, 'Loading events... ');
load(fullfile(fpath_ev, fname_ev));
fprintf(1, 'Done! '); toc;


%% trigger event analyses

% sensor dimensions
width = events.width;
height = events.height;

numPosTriggers = sum(events.trigger.p == 1);
numNegTriggers = sum(events.trigger.p == -1);

ts_FirstTrigger = events.trigger.ts(1);
ts_LastTrigger = events.trigger.ts(end);

posTriggersPerSec = numPosTriggers / ((ts_LastTrigger - ts_FirstTrigger)./1e6);
negTriggersPerSec = numNegTriggers / ((ts_LastTrigger - ts_FirstTrigger)./1e6);

if numPosTriggers ~= length(dd)
    warning('Unequal number of positive triggers and FLIR frames!');
else
    fprintf(1, 'Number of positive triggers equals number of FLIR frames!\n');
end
pause(3);

%% event data subsampling

% add 1 to events due to Matlab indexing
events.x = events.x + 1;
events.y = events.y + 1;

% only consider events within first and last trigger
idx_sub = events.ts >= events.trigger.ts(1) & events.ts <= events.trigger.ts(end);

events_sub.ts = events.ts(idx_sub);
events_sub.x = events.x(idx_sub);
events_sub.y = events.y(idx_sub);
events_sub.p = events.p(idx_sub);


flag_IE_filter = true;
if flag_IE_filter
    tic;

    % twindow_IE_filter = 50000;
    twindow_IE_filter = 10000;

    fprintf(1, ['Event filtering with ', num2str(twindow_IE_filter), ' microsecond window... ']);

    [IE, IEm] = IE_filter(events, twindow_IE_filter);

%     IE = IE & (IEm > 1);

    events.x = events.x(IE);
    events.y = events.y(IE);
    events.p = events.p(IE);
    events.ts = events.ts(IE);
    events.IEm = IEm(IE);

    % figure;
    % histogram(IEm(IE), 'Normalization', 'probability');
    % xlabel('Inceptive Event Magnitude (TE)');
    % ylabel('Probability');


    fprintf(1,'Done! '); toc;
else
    twindow_IE_filter = [];
end



clear('idx_sub');

% set subsampled set of events to have time
events_sub.ts = events_sub.ts - events_sub.ts(1);

% use another subsample of events (optional)
% find positive and negative events for plotting
idx_sub = ones(size(events_sub.ts));
% idx_sub = events.ts <= 4e6;
idx_sub_pos = idx_sub & events_sub.p == 1;
idx_sub_neg = idx_sub & events_sub.p == -1;

% figure;
% scatter3(events_sub.x(idx_sub_pos), events_sub.y(idx_sub_pos), events_sub.ts(idx_sub_pos)./1e6, 1, 'red');
% hold on;
% scatter3(events_sub.x(idx_sub_neg), events_sub.y(idx_sub_neg), events_sub.ts(idx_sub_neg)./1e6, 1, 'blue');
% % zlim([0, 6]);
% view(90, 0);
% title('Positive Event: Red; Negative Event: Blue');
% xlabel('x');
% ylabel('y');
% % zlabel('time (\mus)');
% zlabel('time (s)');
% hold off;

%%


sig_max = squeeze(max(imgs, [], [1,2]));

figure;
plot(tvec_FLIR, sig_max, 'LineWidth', 2.5);
xlim([0, 6]);
% ylim([0, 255]);

%%
% 
% ts_trigger = (events.trigger.ts - events.trigger.ts(1))./1e6;
% idx_trigger = events.trigger.p == 1;
% 
% % figure;
% % plot(tvec_FLIR, sig_max, 'LineWidth', 2.5);
% % hold on;
% % xline(ts_trigger(idx_trigger), 'k--');
% % % xline((events.trigger.ts(1:100) - events.trigger.ts(1))./1e6, 'k--');
% % % xline((events.trigger.ts - events.trigger.ts(1))./1e6, 'k--', 'LineWidth', 0.5);
% % xlim([0, 3]);
% % % ylim([0, 255]);
% % hold off;
% 
% 
% tint = (1/evFrames_FPS)*1e6;
% twin = (1/evFrames_FPS)*1e6; % (+/- twin/2)
% tvec_ev = events_sub.ts(1) : tint : events_sub.ts(end);
% 
% % pre-allocate event frames
% B_arr = zeros(height, width, length(tvec_ev));
% 
% for ii = 1:length(tvec_ev)
%     fprintf(1, ['Generating event frame ', num2str(ii), ' of ', num2str(length(tvec_ev)), '.\n']);
%     twin_tmp = [tvec_ev(ii) - twin/2, tvec_ev(ii) + twin/2];
%     idx = (events_sub.ts >= twin_tmp(1)) & (events_sub.ts < twin_tmp(2));
%     B = accumarray([events_sub.y(idx), events_sub.x(idx)], events_sub.p(idx), [height, width], @sum, 0);
% 
%     B_arr(:,:,ii) = B;
% 
% 
% end
% 
% clear('twin_tmp', 'B', 'idx');
% 
% numEvents = squeeze(sum(B_arr, [1,2]));
% figure;
% plot(tvec_ev./1e6, numEvents./max(numEvents(:)));


%%
n_upsample = evFrames_FPS / FLIR_FPS; % upsample factor
tvec_ev = events.trigger.ts(events.trigger.p==1); % set event data time vector to positive trigger events

firstframe_keep = 1;
tvec_ev = tvec_ev(firstframe_keep:end);

tvec_ev = upsample(tvec_ev, n_upsample); % upsample time vector
tvec_ev = tvec_ev(1:end-n_upsample+1); % remove trailing zeros
idx_tvec_ev = find(tvec_ev); % non-zero indices
idx_tvec_interp = find(tvec_ev == 0); % indices for interpolation of FLIR data
tvec_ev = interp1(idx_tvec_ev, tvec_ev(idx_tvec_ev), 1:length(tvec_ev), "linear"); % interpolate times of upsampled indices
% tvec_ev = (tvec_ev - tvec_ev(1)) ./ 1e6; % start events at 0

% tvec_ev = tvec_ev + 5000;

% height and width of event sensor
height = events.height;
width = events.width;


T = zeros(height, width);
T2 = zeros(height, width);
k = 5.0e-5; % larger k: signal decays faster
% k = 30e-5;

% start of data window
time_start = min(events.ts);

% pre-allocate decay frame
evFrames = zeros(height, width, length(tvec_ev));

for ii = 1:length(tvec_ev)
    fprintf(1, ['Generating exponential decay frame ', num2str(ii), ' of ', num2str(length(tvec_ev)), '.\n']);
    
    % time_end = time_start + tint;
    time_end = tvec_ev(ii); % time of the frame
    
    % find events within the frame
    idx = find((events.ts >= time_start) & (events.ts < time_end));
    
    % calculate the decay of each event from the frame start time
    delta = events.p(idx) .* exp(-k * (time_end - events.ts(idx))); % with polarity
%     delta = abs(events.p(idx)) .* exp(-k * (time_end - events.ts(idx))); % without polarity

    % accumulate the events on a per-pixel basis
    B_tmp = accumarray([events.y(idx), events.x(idx)], delta, [height, width], @sum, 0);

    % decay the previous frame and add the events from the new frame
    T = T .* exp(-k * (time_end - time_start)) + B_tmp;

    evFrames(:,:,ii) = T;

    time_start = time_end; % set window start time to the frame just processed
end

clear('B_tmp', 'B_tmp2');

figure; imagesc(evFrames(:,:,1)); axis image; colormap(gray); colorbar; clim([-1, 1]);


numEvents = squeeze(sum(evFrames, [1,2]));

tvec_ev_plt = (tvec_ev - tvec_ev(1)) ./ 1e6;

figure;
plot(tvec_ev_plt, numEvents./max(numEvents(:)));



%%

% tdiff = 0 - 0.01;
% tdiff = 2.58 - 1.34;

tdiff = 0.01;

figure;
plot(tvec_FLIR, sig_max./255.*10 - 0.5, 'LineWidth', 2.5);
hold on;
plot(tvec_ev_plt, numEvents./max(numEvents(:)), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Signal');
legend('Conventional (Raw)', 'Event', 'location', 'northeast');
hold off;


figure;
plot(tvec_FLIR, sig_max./255.*10 - 0.5, 'LineWidth', 2.5);
hold on;
plot(tvec_FLIR + tdiff, sig_max./255.*10 - 0.5, 'LineWidth', 2.5);
plot(tvec_ev_plt, numEvents./max(numEvents(:)), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Signal');
legend('Conventional (Raw)', 'Conventional with Offset', 'Event', 'location', 'northeast');
hold off;


%%

ts_trigger = events.trigger.ts(events.trigger.p==1);
ts_trigger = (ts_trigger - ts_trigger(1)) ./ 1e6;

% figure;
% plot(tvec_FLIR, sig_max./255.*10 - 0.5, 'LineWidth', 2.5);
% hold on;
% plot(tvec_ev_plt, numEvents./max(numEvents(:)), 'LineWidth', 1.5);
% yline(ts_trigger, 'k--');
% xlabel('Time (s)');
% ylabel('Signal');
% legend('Conventional (Raw)', 'Event', 'location', 'northeast');
% hold off;

ev_signal_diff = diff(numEvents./max(numEvents(:)));
idx_ev_tmp = find(ev_signal_diff >= 0.05);

figure;
plot(tvec_ev_plt(idx_ev_tmp-1), ev_signal_diff(idx_ev_tmp), 'o');


figure;
plot(tvec_FLIR, sig_max./255.*10 - 0.5, 'LineWidth', 2.5);
hold on;
plot(tvec_ev_plt, numEvents./max(numEvents(:)), 'LineWidth', 1.5);
% yline(ts_trigger, 'k--');
% xline(tvec_ev_plt(idx_ev_tmp(1:10)), 'r--');
% xline(5, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Signal');
legend('Conventional (Raw)', 'Event', 'location', 'northeast');
hold off;

%%

x_FLIR = [1.10, 2.38, 3.57, 4.73, 5.89, 7.06]; % last frame before signal increase
% x_ev = [1.15275, 2.38955, 3.58033, 4.73708, 5.90183, 7.07259];
x_ev = [1.15075, 2.38755, 3.57832, 4.73508, 5.89783, 7.07059]; % last event frame before signal increase

x_FLIR_down = [1.75, 2.95, 4.15, 5.29, 6.44, 7.75];
x_ev_down = [1.80317, 2.95392, 4.1527, 5.30345, 6.44419, 7.76104];

dt = x_ev - x_FLIR;
dt_down = x_ev_down - x_FLIR_down;

figure;
plot(x_ev, dt, '-o', 'LineWidth', 2);
hold on;
plot(x_ev_down, dt_down, '-o', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Offset (s)');
legend('Positive Trigger', 'Negative Trigger');
hold off;
