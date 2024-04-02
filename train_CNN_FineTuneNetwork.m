clear; close all; clc;

rng_seed = rng("default");

%% define network


dH5_net = load('./Trained_Networks\LargeNetwork_SlopeMSE_CartesianMeshgrid\depth4\trained_network_depth4_slopeMSE_cartesian.mat');

net = dH5_net.net;

% analyzeNetwork(net);

%% define training parameters

% learnRate_arr = [0.005, 0.005 * 0.25, 0.005 * 0.25 * 0.25];
% learnRate_arr = 0.005 * 0.25 * 0.25;
learnRate_arr = 0.00001;

% learnRate = 0.001; % global learning rate (default = 0.001)
gradDecay = 0.9; % gradient decay factor (default = 0.9)
sqGradDecay = 0.999; % squared gradient decay factor (default = 0.999)
epsilon = 1e-8; % small constant preventing divide-by-zero (default = 1e-8)

l2Regularization = 0.0001; % weight regularization

% max number of epochs
maxEpochs = 1;

% use least-squares lost function starting on this epoch
epoch_use_customCostFunction = 2;

flag_useAdam = true;
if flag_useAdam
    % parameters for Adam optimizer
    % averageGrad = []; % start from scratch
    % averageSqGrad = []; % start from scratch
    averageGrad = dH5_net.averageGrad; % start from pre-trained network
    averageSqGrad = dH5_net.averageSqGrad; % start from pre-trained network

    fprintf(1, ['Training ', num2str(maxEpochs), ' epochs with Adam optimizer.\n']);

else
    % parameters for stochastic gradient descent with momentum (SGDM)
    vel = [];
    momentum = 0.9;

    fprintf(1, ['Training ', num2str(maxEpochs), ' epochs with SGDM optimizer.\n']);

end

if epoch_use_customCostFunction > maxEpochs
    fprintf(1, 'Only using slope MSE loss function.\n');
else
    fprintf(1, ['Using least-squares loss function starting on epoch ', num2str(epoch_use_customCostFunction), '.\n']);
end


%% set up minibatches for Train data

filename_arr_Train = dH5_net.filename_arr_Train;
filename_arr_Valid = dH5_net.filename_arr_Valid;
idx_current_minibatch_Train = dH5_net.idx_current_minibatch_Train;
idx_current_minibatch_Valid = dH5_net.idx_current_minibatch_Valid;
maxIterations = maxEpochs * length(idx_current_minibatch_Train);

%% Set up weighting matrix for least-squares loss function

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

    if epoch >= epoch_use_customCostFunction
        flag_loss_function = 1; % use least-squares
        loss_fcn_str = "Least-Squares";

        % flag_loss_function = 2; % use Strehl
        % loss_fcn_str = "Strehl";
    else
        flag_loss_function = 0; % use slope MSE
        loss_fcn_str = "Slope MSE";
    end

    %%%%%%%%%%%%%%% PUT TRAINING CODE HERE %%%%%%%%%%%%%%%%%%%

    for jj = 1:length(idx_current_minibatch_Train)-1 % loop through all data
    % for jj = 1:73 % specific number of iterations

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

    % loss_Valid_arr = [];
    % 
    % % loop over the datasets
    % for kk = 1:length(idx_current_minibatch_Valid)-1
    % 
    %     % check if "stop" button has been selected
    %     if monitor.Stop
    %         break
    %     end
    % 
    %     idx_load = idx_current_minibatch_Valid(kk):idx_current_minibatch_Valid(kk+1)-1;
    % 
    %     fds_Valid = fileDatastore(filename_arr_Valid(idx_load), 'ReadFcn', ...
    %         @readData_Valid_Test, 'IncludeSubfolders', true, 'FileExtensions', '.mat');
    % 
    %     mbq_Valid = minibatchqueue(fds_Valid, ...
    %         MiniBatchSize=length(fds_Valid.Files), ...
    %         MiniBatchFormat={'SSCB', 'CB', 'B', 'B', 'B'});
    % 
    %     [X, T, dataset_flag, subaperture_flag, frame_flag] = next(mbq_Valid);
    % 
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    %     if flag_loss_function == 0
    % 
    %         W_process = [];
    % 
    %     elseif flag_loss_function == 1 % want to use least-squares
    % 
    %         clear('i', 'idx_tmp', 'W_process', 'dataset_process');
    % 
    %         dataset_tmp = gather(extractdata(dataset_flag));  % all datasets loaded
    %         subaperture_tmp2 = gather(extractdata(subaperture_flag));  % all subapertures loaded
    % 
    %         idx_new_dataset = find(subaperture_tmp2 == 1);  % find where new dataset begins
    %         dataset_tmp2 = dataset_tmp(idx_new_dataset);
    % 
    %         for i = 1:length(dataset_tmp2)
    %             idx_tmp = find(dataset_tmp2(i) == dataset_flag_Reconstructor);
    %             W_process{i, 1} = W{idx_tmp, 1};
    %             dataset_process(i, 1) = dataset_flag_Reconstructor(idx_tmp);
    %         end
    % 
    %         % loss_tmp = modelLoss_Valid_Test(net, X, T, W_process);
    % 
    %     end
    % 
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    %     loss_tmp = modelLoss_Valid_Test(net, X, T, flag_loss_function, W_process);
    % 
    %     if isnan(loss_tmp)
    %         disp("NaN!");
    %     end
    % 
    %     loss_Valid_arr = [loss_Valid_arr; loss_tmp];
    % 
    % end
    % 
    % loss_Valid = mean(loss_Valid_arr);
    % 
    % recordMetrics(monitor, iteration, ValidationLoss=loss_Valid);
    % updateInfo(monitor, ValidationLoss=loss_Valid);

end 

tEnd = toc(tStart);

monitor.Status = "Finished";

%%

tic; fprintf(1, 'Saving... ');
save('trained_network_depth4_slopeMSE_cartesian_FineTuneSlope_test.mat');
fprintf(1, 'Done. '); toc;

%% functions
function data = readData_Train(filename)

    load(filename, 'X', 'Y', 'dataset_flag', 'subaperture_flag', 'frame_flag');

    % idx_keep = [1, 5, 9, 10]; % depth of 1
    % idx_keep = [1, 2, 5, 6, 9, 10]; % depth of 2
    % idx_keep = [1, 2, 3, 5, 6, 7, 9, 10]; % depth of 3
    % idx_keep = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]; % (max) depth of 4
    % 
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    idx_keep = [1, 2, 3, 5, 6, 7, 9, 10]; % depth of 3
    % idx_keep = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]; % (max) depth of 4
    % 
    X = X(:, :, idx_keep);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    % M1 = -30e-3 / 1500e-3; % front 4-f system
    % M2 = -100e-3 / 40e-3; % relay lens
    % f_MLA = 24e-3; % focal length of MLA
    % pixel_pitch_ev = 15e-6; % 15 micrometer pitch event pixels
    % pix2rad = pixel_pitch_ev .* abs(M1) ./ (abs(M2) .* f_MLA);
    pix2rad = 5e-6; % avoid computing every time

    wvl = 635e-9; % wavelength for Strehl computation

    phi_Y = [];
    phi_T = [];

    idx = 1;
    for ii = 1:length(W)
        W2 = dlarray(W{ii}, "");  % least-squares weighting for current dataset
        idx2 = idx:idx+size(W{ii},2)/2-1;

        T2 = [T(1, idx2), T(2, idx2)];
        Y2 = [Y(1, idx2), Y(2, idx2)];
        T2 = (T2.stripdims)'; Y2 = (Y2.stripdims)';

        T2 = T2 .* pix2rad;
        Y2 = Y2 .* pix2rad;

        phi_Y = [phi_Y; reconstruct_phase(W2, Y2, wvl)];
        phi_T = [phi_T; reconstruct_phase(W2, T2, wvl)];

        % phi_Y = reconstruct_phase(W2, Y2, wvl);
        % phi_T = reconstruct_phase(W2, T2, wvl);
        % 
        % [S(ii), ~] = Strehl_approximation(phi_Y, phi_T);

        idx = idx + size(W{ii},2)/2;  % increment index
    end

    S = -mse(phi_T, phi_Y, "DataFormat", "B");
    S = S.exp;
    loss = 1-S;

    % loss = mean(1-S);  % minus sign because minimizing 1-Strehl but derivative of 1 does not matter

end

% Calculate gradients of loss with respect to learnable parameters.
gradients = dlgradient(loss, net.Learnables);
% figure; imagesc(extractdata(gradients.Value{1}(:,:,1,1))); axis image; colorbar; title('Gradients');

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

elseif flag_loss_function == 2  % Strehl

    % M1 = -30e-3 / 1500e-3; % front 4-f system
    % M2 = -100e-3 / 40e-3; % relay lens
    % f_MLA = 24e-3; % focal length of MLA
    % pixel_pitch_ev = 15e-6; % 15 micrometer pitch event pixels
    % pix2rad = pixel_pitch_ev .* abs(M1) ./ (abs(M2) .* f_MLA);
    pix2rad = 5e-6; % avoid computing every time

    wvl = 635e-9; % wavelength for Strehl computation

    phi_Y = [];
    phi_T = [];

    idx = 1;
    for ii = 1:length(W)
        W2 = dlarray(W{ii}, "");  % least-squares weighting for current dataset
        idx2 = idx:idx+size(W{ii},2)/2-1;

        T2 = [T(1, idx2), T(2, idx2)];
        Y2 = [Y(1, idx2), Y(2, idx2)];
        T2 = (T2.stripdims)'; Y2 = (Y2.stripdims)';

        T2 = T2 .* pix2rad;
        Y2 = Y2 .* pix2rad;

        phi_Y = [phi_Y; reconstruct_phase(W2, Y2, wvl)];
        phi_T = [phi_T; reconstruct_phase(W2, T2, wvl)];

        % phi_Y = reconstruct_phase(W2, Y2, wvl);
        % phi_T = reconstruct_phase(W2, T2, wvl);
        % 
        % [S(ii), ~] = Strehl_approximation(phi_Y, phi_T);

        idx = idx + size(W{ii},2)/2;  % increment index
    end

    S = -mse(phi_T, phi_Y, "DataFormat", "B");
    S = S.exp;
    loss = 1-S;

    % loss = mean(1-S);  % minus sign because minimizing 1-Strehl but derivative of 1 does not matter

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