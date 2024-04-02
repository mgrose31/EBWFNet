clear; close all; clc;

rng_seed = rng("default");

%% use for implementing functionality for loading

% fds_Train = fileDatastore('./Datastore/Train/', 'ReadFcn', ...
%     @readData_Train, 'IncludeSubfolders', true, ...
%     'FileExtensions', '.mat');
% fds_Valid = fileDatastore('./Datastore/Validation/', 'ReadFcn', ...
%     @readData_Train, 'IncludeSubfolders', true, ...
%     'FileExtensions', '.mat');
% fds_Test = fileDatastore('./Datastore/Test/', 'ReadFcn', ...
%     @readData_Train, 'IncludeSubfolders', true, ...
%     'FileExtensions', '.mat');
% 
% data = readData_Train(fds_Train.Files{1});
% 
% X_isnan = zeros(size(fds_Train.Files));
% Y_isnan = zeros(size(fds_Train.Files));
% for ii = 1:length(fds_Train.Files)
%     if rem(ii, 100) == 0
%         disp(ii);
%     end
%     data = readData_Train(fds_Train.Files{ii});
%     X_isnan(ii) = any(isnan(data{1}(:)));
%     Y_isnan(ii) = any(isnan(data{2}(:)));
% end


%% define network

% good network
layers = [
    % imageInputLayer([50, 50, 4], 'Normalization', 'none') % depth of 1
    % imageInputLayer([50, 50, 6], 'Normalization', 'none') % depth of 2
    % imageInputLayer([50, 50, 8], 'Normalization', 'none') % depth of 3
    imageInputLayer([50, 50, 10], 'Normalization', 'none') % (max) depth of 4

    convolution2dLayer(15, 10)
    batchNormalizationLayer
    reluLayer

    convolution2dLayer(9, 20)  % default network
    % convolution2dLayer(9, 15)  % small network
    batchNormalizationLayer
    reluLayer

    convolution2dLayer(5, 30)  % default network
    % convolution2dLayer(5, 20)  % small network
    batchNormalizationLayer
    reluLayer

    maxPooling2dLayer(2, 'Stride', 2)

    convolution2dLayer(3, 40)  % default network
    % convolution2dLayer(3, 25)  % small network network
    batchNormalizationLayer
    reluLayer

    maxPooling2dLayer(2, 'Stride', 2)

    convolution2dLayer(3, 50)  % default network
    % convolution2dLayer(3, 30)  % small network
    batchNormalizationLayer
    reluLayer

    maxPooling2dLayer(3, 'Stride', 1)

    fullyConnectedLayer(2)

    % regressionLayer
    ];

net = dlnetwork(layers);

analyzeNetwork(net);

%% define training parameters

learnRate_arr = [0.005, 0.005 * 0.25, 0.005 * 0.25 * 0.25];

% learnRate = 0.001; % global learning rate (default = 0.001)
gradDecay = 0.9; % gradient decay factor (default = 0.9)
sqGradDecay = 0.999; % squared gradient decay factor (default = 0.999)
epsilon = 1e-8; % small constant preventing divide-by-zero (default = 1e-8)

l2Regularization = 0.0001; % weight regularization

% max number of epochs
maxEpochs = 3;

% use WLS / Strehl cost function starting on this epoch
epoch_use_customLossFunction = 4;

flag_useAdam = true;
if flag_useAdam
    % parameters for Adam optimizer
    averageGrad = [];
    averageSqGrad = [];

    fprintf(1, ['Training ', num2str(maxEpochs), ' epochs with Adam optimizer.\n']);

else
    % parameters for stochastic gradient descent with momentum (SGDM)
    vel = [];
    momentum = 0.9;

    fprintf(1, ['Training ', num2str(maxEpochs), ' epochs with SGDM optimizer.\n']);

end

if epoch_use_customLossFunction > maxEpochs
    fprintf(1, 'Only using slope MSE loss function.\n');
else
    fprintf(1, ['Using least-squares loss function starting on epoch ', num2str(epoch_use_customLossFunction), '.\n']);
end


%% set up minibatches for Train data

% dd_Train = dir('./Datastore/Train/dataset*');  % 5 ms IE
% dd_Valid = dir('./Datastore/Validation/dataset*');  % 5 ms IE
% dd_Train = dir('C:/Users/ISSL/Documents/Datastore_CNN_2p5msIE/Train/dataset*');  % 2.5 ms IE
% dd_Valid = dir('C:/Users/ISSL/Documents/Datastore_CNN_2p5msIE/Validation/dataset*');  % 2.5 ms IE
dd_Train = dir('C:/Users/ISSL/Documents/Datastore_CNN_10msIE/Train/dataset*');  % 10 ms IE
dd_Valid = dir('C:/Users/ISSL/Documents/Datastore_CNN_10msIE/Validation/dataset*');  % 10 ms IE
% dd_Train = dir('./Datastore_noIE/Train/dataset*');  % no IE
% dd_Valid = dir('./Datastore_noIE/Validation/dataset*');  % no IE

dataset_arr_Train = [];
frame_arr_Train = [];
dd_arr_Train = [];

for ii = 1:length(dd_Train)
% for ii = 1
    dataset_flags_Train_tmp(ii, 1) = str2double(dd_Train(ii).name(end-1:end));
    dd1 = dir(fullfile(dd_Train(ii).folder, dd_Train(ii).name, '/Frame*'));
    numFrames_Train(ii, 1) = length(dd1);
    dataset_arr_Train = [dataset_arr_Train; ones(numFrames_Train(ii), 1) * dataset_flags_Train_tmp(ii)];
    frame_arr_Train = [frame_arr_Train; (1:numFrames_Train(ii))'];
    dd_arr_Train = [dd_arr_Train; dd1];
end
clear('dd1');

for ii = 1:length(dd_arr_Train)
    filename_arr_Train{ii, 1} = fullfile(dd_arr_Train(ii).folder, dd_arr_Train(ii).name);
end

idx_randomized_dataset_frame = randperm(length(filename_arr_Train));

filename_arr_Train = filename_arr_Train(idx_randomized_dataset_frame);
dataset_arr_Train = dataset_arr_Train(idx_randomized_dataset_frame);
frame_arr_Train = frame_arr_Train(idx_randomized_dataset_frame);

idx_current_minibatch_Train = 1:8:length(filename_arr_Train);
if idx_current_minibatch_Train(end) ~= length(filename_arr_Train)
    idx_current_minibatch_Train = [idx_current_minibatch_Train, length(filename_arr_Train)];
end

maxIterations = maxEpochs * length(idx_current_minibatch_Train);

%% set up minibatches for Validation data

dataset_arr_Valid = [];
frame_arr_Valid = [];
dd_arr_Valid = [];

for ii = 1:length(dd_Valid)
% for ii = 1
    dataset_flags_Valid_tmp(ii, 1) = str2double(dd_Valid(ii).name(end-1:end));
    dd1 = dir(fullfile(dd_Valid(ii).folder, dd_Valid(ii).name, '/Frame*'));
    numFrames_Valid(ii, 1) = length(dd1);
    dataset_arr_Valid = [dataset_arr_Valid; ones(numFrames_Valid(ii), 1) * dataset_flags_Valid_tmp(ii)];
    frame_arr_Valid = [frame_arr_Valid; (1:numFrames_Valid(ii))'];
    dd_arr_Valid = [dd_arr_Valid; dd1];
end
clear('dd1');

for ii = 1:length(dd_arr_Valid)
    filename_arr_Valid{ii, 1} = fullfile(dd_arr_Valid(ii).folder, dd_arr_Valid(ii).name);
end

idx_current_minibatch_Valid = 1:8:length(filename_arr_Valid);
if idx_current_minibatch_Valid(end) ~= length(filename_arr_Valid)
    idx_current_minibatch_Valid = [idx_current_minibatch_Valid, length(filename_arr_Valid)];
end

%% Set yp weighting matrix for least-squares loss function

dd_Geometry = dir('./WFS_Reconstruction/Reconstructors/*.mat');

for ii = 1:length(dd_Geometry)
    d_Geometry = load(fullfile(dd_Geometry(ii).folder, dd_Geometry(ii).name));

    Reconstructor{ii, 1} = d_Geometry.Reconstructor_TiltIncluded;  % use tilt-included reconstructor
    % Reconstructor{ii, 1} = d_Geometry.Reconstructor_TiltRemoved;  % use tilt-removed reconstructor

    W{ii, 1} = Reconstructor{ii, 1}' * Reconstructor{ii, 1};  % least-squares weighting

    dataset_flag_Reconstructor(ii, 1) = str2double(dd_Geometry(ii).name(end-5:end-4));  % keep track of dataset

    % figure;
    % imagesc(W{ii, 1});
    % axis image;
    % colorbar;

end

%% test getting the correct least-squares weighting

% fds_Train = fileDatastore('./Datastore/Train/', ...
%     'ReadFcn', @myReadFcn_Train, 'IncludeSubfolders', true, ...
%     'FileExtensions', '.mat');
% % data_tmp = read(fds_Train);
% 
% mbq = minibatchqueue(fds_Train, ...
%     MiniBatchSize=8, ...
%     MiniBatchFormat={'SSCB', 'CB', 'B', 'B', 'B'});
% 
% shuffle(mbq);


% idx_load = idx_current_minibatch(1):idx_current_minibatch(1+1)-1;
% 
% fds_Train = fileDatastore(filename_arr_Train(idx_load), 'ReadFcn', ...
%     @myReadFcn_Train, 'IncludeSubfolders', true, 'FileExtensions', '.mat');
% 
% mbq = minibatchqueue(fds_Train, ...
%     MiniBatchSize=length(fds_Train.Files), ...
%     MiniBatchFormat={'SSCB', 'CB', 'B', 'B', 'B'});
% 
% [Xtmp, Ttmp, dataset_tmp, subaperture_tmp, frame_tmp] = next(mbq);
% 
% dataset_tmp2 = gather(extractdata(dataset_tmp));  % all datasets loaded
% subaperture_tmp2 = gather(extractdata(subaperture_tmp));  % all subapertures loaded
% 
% idx_new_dataset = find(subaperture_tmp2 == 1);  % find where new dataset begins
% dataset_tmp3 = dataset_tmp2(idx_new_dataset);
% 
% for ii = 1:length(dataset_tmp3)
%     idx_tmp = find(dataset_tmp3(ii) == dataset_flag_Reconstructor);
%     W_process{ii, 1} = W{idx_tmp, 1};
%     dataset_process(ii, 1) = dataset_flag_Reconstructor(idx_tmp);
% end


%% train network

monitor = trainingProgressMonitor( ...
    Metrics = ["TrainingLoss", "ValidationLoss"], ...
    Info = ["Epoch", "LearnRate", "Iteration", "TrainingLoss", "ValidationLoss", "LossFcn"], ...
    XLabel = "Iteration");
groupSubPlot(monitor, "Loss", ["TrainingLoss", "ValidationLoss"]);

monitor.Status = "Training";

tStart = tic;

epoch = 0;
iteration = 0;

while epoch < maxEpochs && ~monitor.Stop

    epoch = epoch + 1;

    learnRate = learnRate_arr(epoch);

    if epoch >= epoch_use_customLossFunction
        flag_loss_function = 1; % use WLS
        loss_fcn_str = "Weighted Least-Squares";

        % flag_loss_function = 2; % use approximate Strehl
        % loss_fcn_str = "Strehl";
    else
        flag_loss_function = 0; % use slope MSE
        loss_fcn_str = "Slope MSE";
    end

    %%%%%%%%%%%%%%% PUT TRAINING CODE HERE %%%%%%%%%%%%%%%%%%%

    for jj = 1:length(idx_current_minibatch_Train)-1

        % check if "stop" button has been selected
        if monitor.Stop
            break
        end

        iteration = iteration + 1;

        monitor.Progress = 100 * iteration / maxIterations;

        idx_load = idx_current_minibatch_Train(jj):idx_current_minibatch_Train(jj+1)-1;

        fds_Train = fileDatastore(filename_arr_Train(idx_load), 'ReadFcn', ...
            @readData_Train, 'IncludeSubfolders', true, 'FileExtensions', '.mat');

        mbq = minibatchqueue(fds_Train, ...
            MiniBatchSize=length(fds_Train.Files), ...
            MiniBatchFormat={'SSCB', 'CB', 'B', 'B', 'B'});

        [X, T, dataset_flag, subaperture_flag, frame_flag] = next(mbq);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if flag_loss_function == 0

            W_process = [];

        elseif flag_loss_function == 1 || flag_loss_function == 2 % weighted least-squares or Strehl

            clear('i', 'idx_tmp', 'W_process', 'dataset_process');

            dataset_tmp = gather(extractdata(dataset_flag));  % all datasets loaded
            subaperture_tmp2 = gather(extractdata(subaperture_flag));  % all subapertures loaded

            idx_new_dataset = find(subaperture_tmp2 == 1);  % find where new dataset begins
            dataset_tmp2 = dataset_tmp(idx_new_dataset);

            if flag_loss_function == 1 % weighted least-squares

                for i = 1:length(dataset_tmp2)
                    idx_tmp = find(dataset_tmp2(i) == dataset_flag_Reconstructor);
                    W_process{i, 1} = W{idx_tmp, 1}; % grab W
                    dataset_process(i, 1) = dataset_flag_Reconstructor(idx_tmp);
                end

            elseif flag_loss_function == 2 % Strehl

                for i = 1:length(dataset_tmp2)
                    idx_tmp = find(dataset_tmp2(i) == dataset_flag_Reconstructor);
                    W_process{i, 1} = Reconstructor{idx_tmp, 1}; % grab reconstructor
                    dataset_process(i, 1) = dataset_flag_Reconstructor(idx_tmp);
                end

            else

                error("Error!");

            end
            % loss_tmp = modelLoss_Valid_Test(net, X, T, flag_loss_function, W_process);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [loss, gradients, state] = dlfeval(@modelLoss_Train, net, X, T, flag_loss_function, W_process);
        net.State = state;

        if isnan(loss)
            disp("NaN!");
        end

        % apply L2 Regularization to the weights
        idx = net.Learnables.Parameter == "Weights";
        gradients(idx,:) = dlupdate(@(g,w) g + l2Regularization*w, gradients(idx,:), net.Learnables(idx,:));

        if flag_useAdam
            [net, averageGrad, averageSqGrad] = adamupdate(net, gradients, ...
                averageGrad, averageSqGrad, iteration, ...
                learnRate, gradDecay, sqGradDecay, epsilon);
        else
            [net, vel] = sgdmupdate(net, gradients, vel, learnRate, momentum);
        end

        recordMetrics(monitor, iteration, TrainingLoss=loss);
        updateInfo(monitor, ...
            Epoch=string(epoch) + " of " + string(maxEpochs), ...
            Iteration=iteration, ...
            LearnRate=learnRate, ...
            TrainingLoss=string(loss), ...
            LossFcn=loss_fcn_str);

    end

    %%%%%%% VALIDATION CODE HERE %%%%%%%%

    loss_Valid_arr = [];

    % loop over the datasets
    for kk = 1:length(idx_current_minibatch_Valid)-1

        % check if "stop" button has been selected
        if monitor.Stop
            break
        end

        idx_load = idx_current_minibatch_Valid(kk):idx_current_minibatch_Valid(kk+1)-1;

        fds_Valid = fileDatastore(filename_arr_Valid(idx_load), 'ReadFcn', ...
            @readData_Valid_Test, 'IncludeSubfolders', true, 'FileExtensions', '.mat');

        mbq_Valid = minibatchqueue(fds_Valid, ...
            MiniBatchSize=length(fds_Valid.Files), ...
            MiniBatchFormat={'SSCB', 'CB', 'B', 'B', 'B'});

        [X, T, dataset_flag, subaperture_flag, frame_flag] = next(mbq_Valid);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if flag_loss_function == 0

            W_process = [];

        elseif flag_loss_function == 1 % want to use least-squares

            clear('i', 'idx_tmp', 'W_process', 'dataset_process');

            dataset_tmp = gather(extractdata(dataset_flag));  % all datasets loaded
            subaperture_tmp2 = gather(extractdata(subaperture_flag));  % all subapertures loaded

            idx_new_dataset = find(subaperture_tmp2 == 1);  % find where new dataset begins
            dataset_tmp2 = dataset_tmp(idx_new_dataset);

            for i = 1:length(dataset_tmp2)
                idx_tmp = find(dataset_tmp2(i) == dataset_flag_Reconstructor);
                W_process{i, 1} = W{idx_tmp, 1};
                dataset_process(i, 1) = dataset_flag_Reconstructor(idx_tmp);
            end

            % loss_tmp = modelLoss_Valid_Test(net, X, T, W_process);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        loss_tmp = modelLoss_Valid_Test(net, X, T, flag_loss_function, W_process);

        if isnan(loss_tmp)
            disp("NaN!");
        end

        loss_Valid_arr = [loss_Valid_arr; loss_tmp];

    end

    loss_Valid = mean(loss_Valid_arr);

    recordMetrics(monitor, iteration, ValidationLoss=loss_Valid);
    updateInfo(monitor, ValidationLoss=loss_Valid);

end 

tEnd = toc(tStart);

monitor.Status = "Finished";
monitor.Stop;

%%

tic; fprintf(1, 'Saving... ');
save('trained_network_depth4_slopeMSE_cartesian_10msIE.mat');
fprintf(1, 'Done. '); toc;

%% functions
function data = readData_Train(filename)

    load(filename, 'X', 'Y', 'dataset_flag', 'subaperture_flag', 'frame_flag');

    % idx_keep = [1, 5, 9, 10]; % depth of 1
    % idx_keep = [1, 2, 5, 6, 9, 10]; % depth of 2
    % idx_keep = [1, 2, 3, 5, 6, 7, 9, 10]; % depth of 3
    % idx_keep = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]; % (max) depth of 4

    % idx_keep = 1:8; % no mesh grids

    % X = X(:, :, idx_keep);

    % minTime = 150; %any amount less than 150 microseconds can be ignored (helps with log scaling) (feature normalization)
    % maxTime = 5e6; %any amount greater than 5 seconds can be ignored (put data on fixed output size) (feature normalization)
    % fill_value = single(log(maxTime+1) - log(minTime+1));
    
    fill_value = single(log(5e6 + 1) - log(150 + 1)); % fill value for TORE frames
    v = randi([-10, 10], 1, 2); % generate random integer translations

    % X = imtranslate(X, v, 'linear', 'OutputView', 'same', 'FillValues', fill_value); 
    X(:, :, 1:end-2) = imtranslate(X(:, :, 1:end-2), v, 'linear', 'OutputView', 'same', 'FillValues', fill_value); % only translate TORE frames
    % Y2 = Y + v; % for plotting
    Y = Y + v; % apply translation to output

    % figure;
    % imagesc(X(1,:,end-1), X(:,1,end), X(:,:,1));
    % axis image;
    % colormap(flipud(gray)); colorbar; clim([0, fill_value]);
    % hold on;
    % plot(Y(1), Y(2), 'rx', 'LineWidth', 2);
    % 
    % figure;
    % imagesc(X(1,:,end-1), X(:,1,end), X(:,:,5));
    % axis image;
    % colormap(flipud(gray)); colorbar; clim([0, fill_value]);
    % hold on;
    % plot(Y(1), Y(2), 'rx', 'LineWidth', 2);
    % 
    % figure;
    % imagesc(X(1,:,end-1), X(:,1,end), X(:,:,end-1));
    % axis image;
    % colormap(gray); colorbar; %clim([0, fill_value]);
    % hold on;
    % plot(Y(1), Y(2), 'rx', 'LineWidth', 2);
    % 
    % figure;
    % imagesc(X(1,:,end-1), X(:,1,end), X(:,:,end));
    % axis image;
    % colormap(gray); colorbar; %clim([0, fill_value]);
    % hold on;
    % plot(Y(1), Y(2), 'rx', 'LineWidth', 2);
 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % uncomment for polar coordinates
    % 
    % xx = X(:,:,end-1); % xx meshgrid
    % yy = X(:,:,end); % yy meshgrid
    % [theta, r] = cart2pol(xx, yy); % theta & r meshgrids
    % 
    % [theta_T, r_T] = cart2pol(Y(1), Y(2)); % spot center position in theta & r
    % 
    % X(:, :, end-1) = r; % replace xx meshgrid with r
    % X(:, :, end) = theta; % replace yy meshgrid with theta
    % 
    % % figure;
    % % imagesc(xx(1,:), yy(:,1), theta);
    % % axis image xy square;
    % % colorbar;
    % % hold on;
    % % plot(Y(1), Y(2), 'rx', 'LineWidth', 2);
    % % xlabel('x [m]');
    % % ylabel('y [m]');
    % % hold off;
    % % 
    % % figure;
    % % imagesc(xx(1,:), yy(:,1), r);
    % % axis image xy square;
    % % colorbar;
    % % hold on;
    % % plot(Y(1), Y(2), 'rx', 'LineWidth', 2);
    % % xlabel('x [m]');
    % % ylabel('y [m]');
    % % hold off;
    % 
    % Y = [r_T, theta_T]; % output is r, theta

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    data = [{X}, {Y'}, {dataset_flag}, {subaperture_flag}, {frame_flag}];

end

function data = readData_Valid_Test(filename)

    load(filename, 'X', 'Y', 'dataset_flag', 'subaperture_flag', 'frame_flag');

    % idx_keep = [1, 5, 9, 10]; % depth of 1
    % idx_keep = [1, 2, 5, 6, 9, 10]; % depth of 2
    % idx_keep = [1, 2, 3, 5, 6, 7, 9, 10]; % depth of 3
    % idx_keep = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]; % (max) depth of 4

    % idx_keep = 1:8; % no mesh grids

    % X = X(:, :, idx_keep);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % uncomment for polar coordinates
    % 
    % xx = X(:,:,end-1); % xx meshgrid
    % yy = X(:,:,end); % yy meshgrid
    % [theta, r] = cart2pol(xx, yy); % theta & r meshgrids
    % 
    % [theta_T, r_T] = cart2pol(Y(1), Y(2)); % spot center position in theta & r
    % 
    % X(:, :, end-1) = r; % replace xx meshgrid with r
    % X(:, :, end) = theta; % replace yy meshgrid with theta
    % 
    % % figure;
    % % imagesc(xx(1,:), yy(:,1), theta);
    % % axis image xy square;
    % % colorbar;
    % % hold on;
    % % plot(Y(1), Y(2), 'rx', 'LineWidth', 2);
    % % xlabel('x [m]');
    % % ylabel('y [m]');
    % % hold off;
    % %
    % % figure;
    % % imagesc(xx(1,:), yy(:,1), r);
    % % axis image xy square;
    % % colorbar;
    % % hold on;
    % % plot(Y(1), Y(2), 'rx', 'LineWidth', 2);
    % % xlabel('x [m]');
    % % ylabel('y [m]');
    % % hold off;
    % 
    % Y = [r_T, theta_T]; % output is r, theta

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    data = [{X}, {Y}, {dataset_flag}, {subaperture_flag}, {frame_flag}];

end

function [loss, gradients, state] = modelLoss_Train(net, X, T, flag_loss_function, W)

% Forward data through network.
[Y, state] = forward(net, X);

if flag_loss_function == 0 % MSE between target and estimation slopes

    loss = mse(Y, T); % half mean-squared error

elseif flag_loss_function == 1 % least-squares loss function

    idx = 1;
    for ii = 1:length(W)
        W2 = dlarray(W{ii}, "");  % least-squares weighting for current dataset
        idx2 = idx:idx+length(W{ii})/2-1;

        T2 = [T(1, idx2), T(2, idx2)];
        Y2 = [Y(1, idx2), Y(2, idx2)];
        T2 = T2.stripdims; Y2 = Y2.stripdims;

        loss(ii) = (T2 - Y2) * W2 * (T2 - Y2)';  % compute least-squares loss
        idx = idx + length(W{ii})/2;  % increment index
    end

    loss = mean(loss);  % compute average loss

elseif flag_loss_function == 2  % Strehl

    M1 = -30e-3 / 1500e-3; % front 4-f system
    M2 = -100e-3 / 40e-3; % relay lens
    f_MLA = 24e-3; % focal length of MLA

    pixel_pitch_ev = 15e-6; % 15 micrometer pitch event pixels

    wvl = 635e-9; % wavelength for Strehl computation

    pix2rad = pixel_pitch_ev .* abs(M1) ./ (abs(M2) .* f_MLA);


    idx = 1;
    for ii = 1:length(W)
        W2 = dlarray(W{ii}, "");  % Reconstructor
        idx2 = idx:idx+size(W{ii},2)/2-1;

        T2 = [T(1, idx2), T(2, idx2)];
        Y2 = [Y(1, idx2), Y(2, idx2)];
        T2 = (T2.stripdims)'; Y2 = (Y2.stripdims)';

        T2 = T2 .* pix2rad;
        Y2 = Y2 .* pix2rad;

        phi_Y = reconstruct_phase(W2, Y2, wvl);
        phi_T = reconstruct_phase(W2, T2, wvl);

        [S(ii), ~] = Strehl_approximation(phi_Y, phi_T);

        idx = idx + size(W{ii},2)/2;  % increment index
    end

    loss = mean(-S);  % minus sign because minimizing 1-Strehl but derivative of 1 does not matter

end

% Calculate gradients of loss with respect to learnable parameters.
gradients = dlgradient(loss, net.Learnables);

end

function [loss] = modelLoss_Valid_Test(net, X, T, flag_loss_function, W)

Y = predict(net, X);

if flag_loss_function == 0 % MSE between target and estimation slopes

    loss = mse(Y.stripdims, T.stripdims', "DataFormat", "CB"); % half mean-squared error

elseif flag_loss_function == 1

    T = (T.stripdims)';
    Y = (Y.stripdims);

    idx = 1;
    for ii = 1:length(W)
        W2 = dlarray(W{ii}, "");  % least-squares weighting for current dataset
        idx2 = idx:idx+length(W{ii})/2-1;

        T2 = [T(1, idx2), T(2, idx2)];
        Y2 = [Y(1, idx2), Y(2, idx2)];
        T2 = T2.stripdims; Y2 = Y2.stripdims;

        loss(ii) = (T2 - Y2) * W2 * (T2 - Y2)';  % compute least-squares loss
        idx = idx + length(W{ii})/2;  % increment index
    end

    loss = mean(loss);  % compute average loss

end

end


function OPD = reconstruct_phase(R, s, wvl)
% R is reconstruction matrix
% s is slope vector (all x slopes then all y slopes)

OPD_m = R * s;  % OPD in meters
OPD = OPD_m * (2*pi) / wvl; % OPD in phase (radians)

% figure;
% plot(OPD_m .* 1e6, '-o');
% xlabel('OPD Index');
% ylabel('OPD (\mum)');

end


function [S, MSE] = Strehl_approximation(OPD1, OPD2)
% OPD1 & OPD2 are phase in radians

MSE = mean((OPD1 - OPD2).^2);
S = exp(-MSE);
end