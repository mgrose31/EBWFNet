clear; close all; clc;

% Code to load event data (converted from .raw to .dat using Prophesee
% functions) and trigger data (if it exists), then write the events and
% supporting meta data to a .mat file for later processing.
% Author: Mitchell Grose, University of Dayton, 2023

%% function paths
addpath('../Loading/');

%% set up filename of data to load

% filepath to data
fpath = 'D:\Data\Dual_WFS\2023-07-30_Eastwood\DualSync\Event\dataset35';

% name of event data, trigger, and bias files
fname = 'EventData_cd.dat';
fname_trigger = 'EventData_trigger.dat';
fname_bias = 'EventData.bias';

% get full filepath
filename = fullfile(fpath, fname);
filename_trigger = fullfile(fpath, fname_trigger);
filename_bias = fullfile(fpath, fname_bias);

%% sensor dimensions

% EVK1 VGA sensor
width = 640; % sensor width
height = 480; % sensor height

% % EVK2 HD sensor
% width = 1280; % sensor width
% height = 720; % sensor height

%% load data

tic; fprintf(1, 'Loading events... ');

% load events
events = load_cd_events(filename, false, false, Inf);

% load trigger events if file exists
if isfile(filename_trigger)
    fprintf(1, 'Loading trigger events... ');
    events_trigger = load_ext_trigger_data(filename_trigger);
    events.trigger.ts = events_trigger.ts;
    events.trigger.p = events_trigger.p;
    events.trigger.id = events_trigger.id;
else
    fprintf(1, 'No trigger events... ');
end

% read biases if biases exist
if isfile(filename_bias)
    fprintf(1, 'Loading biases... ');
    biases = readtable(filename_bias, 'FileType', 'text'); % read file
    biases = removevars(biases, 'Var2'); % remove variable of percentage signs
    biases = renamevars(biases, ["Var1", "Var3"], ["Value", "Name"]); % change variable names
    biases = movevars(biases, "Value", "After", "Name"); % put variable 'Value' after variable 'Name'

    events.biases = biases;
else
    fprintf(1, 'No biases... ');
end

fprintf(1, 'Done! '); toc;

%% report statistics on trigger events

if isfile(filename_trigger)
    numPosTriggers = sum(events.trigger.p==1);
    numNegTriggers = sum(events.trigger.p==-1);

    triggerStartTime = events.trigger.ts(1) / 1e6;
    triggerEndTime = events.trigger.ts(end) / 1e6;

    PosTriggersPerSecond = numPosTriggers / (triggerEndTime - triggerStartTime);

    fprintf(1, [num2str(numPosTriggers), ' positive trigger events.\n']);
    fprintf(1, [num2str(numNegTriggers), ' negative trigger events.\n']);
    fprintf(1, [num2str(numPosTriggers), ' FLIR frames expected.\n']);

    fprintf(1, [num2str(PosTriggersPerSecond), ' positive triggers per second.\n']);

end

%% save as mat file

events.width = width;
events.height = height;

save([filename(1:end-4), '.mat'], 'events', '-v7.3');
