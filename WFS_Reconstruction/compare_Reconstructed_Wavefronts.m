clear; close all; clc;

% path to reconstruction functions
addpath(genpath('./mfiles/'));

%% load data for comparison

dataset_flag = 35;

% recon_str = 'TimeAvgIncluded_FullApTiltIncluded';
% recon_str = 'TimeAvgIncluded_FullApTiltRemoved';
% recon_str = 'TimeAvgRemoved_FullApTiltIncluded'; % use this for estimating r0
recon_str = 'TimeAvgRemoved_FullApTiltRemoved'; % use this for metrics

d_SlopeMSE = load(['.\Reconstructed_Wavefronts\', recon_str, '\LargeNetwork_SlopeMSE_CartesianMeshgrid_10msIE\depth4\dataset', num2str(dataset_flag, '%02d'), '_Reconstructed_Wavefronts.mat']);
d_FineTuneWLS = load(['.\Reconstructed_Wavefronts\', recon_str, '\LargeNetwork_SlopeMSE_CartesianMeshgrid_FineTuneWLS\depth4\dataset', num2str(dataset_flag, '%02d'), '_Reconstructed_Wavefronts.mat']);
d_FineTuneStrehl = load(['.\Reconstructed_Wavefronts\', recon_str, '\LargeNetwork_SlopeMSE_CartesianMeshgrid_FineTuneStrehl\depth4\dataset', num2str(dataset_flag, '%02d'), '_Reconstructed_Wavefronts.mat']);

d_Time = load(['..\Formatted_Data\formatted_dataset', num2str(dataset_flag, '%02d'), '.mat'], ...
    'events', 'frameID', 'positions', 'positions_baseline', 'tvec_ev', ...
    'xc_ref_ev_rounded', 'yc_ref_ev_rounded');

tvec_ev = d_Time.tvec_ev;
tvec_ev = (tvec_ev - tvec_ev(1)) ./ 1e6;

%% compare slopes

pix2rad = d_SlopeMSE.pix2rad;

slopes_Frame = d_SlopeMSE.s_FLIR_radians;
slopes_AR1 = d_SlopeMSE.s_Baseline_radians;
slopes_CNN = d_SlopeMSE.s_CNN_radians;
slopes_CNN_FineTuneWLS = d_FineTuneWLS.s_CNN_radians;
slopes_CNN_FineTuneStrehl = d_FineTuneStrehl.s_CNN_radians;

fig_xSlopes = figure('Position', [100, 100, 1000, 200]);
set(gca, 'FontWeight', 'bold');
hold on;
plot(tvec_ev, slopes_Frame(7,:) .* 1e6, 'k-', 'LineWidth', 2);
plot(tvec_ev, slopes_AR1(7,:) .* 1e6, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
plot(tvec_ev, slopes_CNN(7,:) .* 1e6, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
xlim([11, 14]);
ylim([-25, 25]);
xlabel('Time (s)');
ylabel('x-axis Slope (\murad)');
legend('Frame', 'AR(1)', 'EBWFNet', 'location', 'best');
hold off;

fig_ySlopes = figure('Position', [100, 100, 1000, 200]);
set(gca, 'FontWeight', 'bold');
hold on;
plot(tvec_ev, slopes_Frame(47,:) .* 1e6, 'k-', 'LineWidth', 2);
plot(tvec_ev, slopes_AR1(47,:) .* 1e6, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
plot(tvec_ev, slopes_CNN(47,:) .* 1e6, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
xlim([11, 14]);
ylim([-25, 25]);
xlabel('Time (s)');
ylabel('y-axis Slope (\murad)');
legend('Frame', 'AR(1)', 'EBWFNet', 'location', 'best');
hold off;

% exportgraphics(fig_xSlopes, 'Figure_xSlopes_Dataset04_Subaperture07.pdf');
% exportgraphics(fig_ySlopes, 'Figure_ySlopes_Dataset04_Subaperture07.pdf');

sx_FLIR_fullApTilt = pix2rad .* d_SlopeMSE.sx_FLIR_fullApTilt;
sx_Baseline_fullApTilt = pix2rad .* d_SlopeMSE.sx_Baseline_fullApTilt;
sx_CNN_SlopeMSE_fullApTilt = pix2rad .* d_SlopeMSE.sx_CNN_fullApTilt;

sy_FLIR_fullApTilt = pix2rad .* d_SlopeMSE.sy_FLIR_fullApTilt;
sy_Baseline_fullApTilt = pix2rad .* d_SlopeMSE.sy_Baseline_fullApTilt;
sy_CNN_SlopeMSE_fullApTilt = pix2rad .* d_SlopeMSE.sy_CNN_fullApTilt;

figure;
plot(sx_FLIR_fullApTilt .* 1e6, 'k-', 'LineWidth', 2);
hold on;
plot(sx_Baseline_fullApTilt .* 1e6, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
plot(sx_CNN_SlopeMSE_fullApTilt .* 1e6, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
hold off;

fig_fullApTilt = figure('Position', [100, 100, 1000, 400]);

ax1 = subplot(2, 1, 1);
set(gca, 'FontWeight', 'bold');
hold on;
plot(tvec_ev, sx_FLIR_fullApTilt .* 1e6, 'k-', 'LineWidth', 2);
plot(tvec_ev, sx_Baseline_fullApTilt .* 1e6, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
plot(tvec_ev, sx_CNN_SlopeMSE_fullApTilt .* 1e6, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
% xlim([11, 14]);
xlabel('Time (s)');
ylabel('x-axis Full-Ap. Tilt (\murad)');
legend('FLIR', 'Avg. Position', 'CNN (Ours)', 'location', 'best');
hold off;

ax2 = subplot(2, 1, 2);
set(gca, 'FontWeight', 'bold');
hold on;
plot(tvec_ev, sy_FLIR_fullApTilt .* 1e6, 'k-', 'LineWidth', 2);
plot(tvec_ev, sy_Baseline_fullApTilt .* 1e6, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
plot(tvec_ev, sy_CNN_SlopeMSE_fullApTilt .* 1e6, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
% xlim([11, 14]);
xlabel('Time (s)');
ylabel('y-axis Full-Ap. Tilt (\murad)');
legend('FLIR', 'Avg. Position', 'CNN (Ours)', 'location', 'southeast');
hold off;

linkaxes([ax1, ax2], 'xy');

% exportgraphics(fig_fullApTilt, 'Figure_FullApTilt_Dataset04.pdf');


numSubaps = size(slopes_Frame, 1) / 2;

xSlopes_Frame = slopes_Frame(1:numSubaps,:);
ySlopes_Frame = slopes_Frame(numSubaps+1:end,:);
xSlopes_AR1 = slopes_AR1(1:numSubaps,:);
ySlopes_AR1 = slopes_AR1(numSubaps+1:end,:);
xSlopes_CNN = slopes_CNN(1:numSubaps,:);
ySlopes_CNN = slopes_CNN(numSubaps+1:end,:);
xSlopes_CNN_FineTuneWLS = slopes_CNN_FineTuneWLS(1:numSubaps,:);
ySlopes_CNN_FineTuneWLS = slopes_CNN_FineTuneWLS(numSubaps+1:end,:);
xSlopes_CNN_FineTuneStrehl = slopes_CNN_FineTuneStrehl(1:numSubaps,:);
ySlopes_CNN_FineTuneStrehl = slopes_CNN_FineTuneStrehl(numSubaps+1:end,:);

RSE_Slopes_AR1 = sqrt((xSlopes_AR1 - xSlopes_Frame).^2 + (ySlopes_AR1 - ySlopes_Frame).^2);
RSE_Slopes_CNN = sqrt((xSlopes_CNN - xSlopes_Frame).^2 + (ySlopes_CNN - ySlopes_Frame).^2);
RSE_Slopes_CNN_FineTuneWLS = sqrt((xSlopes_CNN_FineTuneWLS - xSlopes_Frame).^2 + (ySlopes_CNN_FineTuneWLS - ySlopes_Frame).^2);
RSE_Slopes_CNN_FineTuneStrehl = sqrt((xSlopes_CNN_FineTuneStrehl - xSlopes_Frame).^2 + (ySlopes_CNN_FineTuneStrehl - ySlopes_Frame).^2);

% RSE_Slope_Baseline_Median = median(RSE_Slope_Baseline(:));
% RSE_Slope_CNN_SlopeMSE_Median = median(RSE_Slope_CNN_SlopeMSE(:));
% RSE_Slope_CNN_FineTuneWLS_Median = median(RSE_Slope_CNN_FineTuneWLS(:));
% RSE_Slope_CNN_FineTuneStrehl_Median = median(RSE_Slope_CNN_FineTuneStrehl(:));

RSE_Slopes_AR1_Mean = mean(RSE_Slopes_AR1(:));
RSE_Slopes_CNN_Mean = mean(RSE_Slopes_CNN(:)); % expect this to be the smallest (best)
RSE_Slopes_CNN_FineTuneWLS_Mean = mean(RSE_Slopes_CNN_FineTuneWLS(:));
RSE_Slopes_CNN_FineTuneStrehl_Mean = mean(RSE_Slopes_CNN_FineTuneStrehl(:));

fprintf(1, ['Slope RSE Mean, AR1: ', num2str(round(RSE_Slopes_AR1_Mean .* 1e6, 3)), '\n']);
fprintf(1, ['Slope RSE Mean, CNN: ', num2str(round(RSE_Slopes_CNN_Mean .* 1e6, 3)), '\n']);
fprintf(1, ['Slope RSE Mean, WLS: ', num2str(round(RSE_Slopes_CNN_FineTuneWLS_Mean .* 1e6, 3)), '\n']);
fprintf(1, ['Slope RSE Mean, Strehl: ', num2str(round(RSE_Slopes_CNN_FineTuneStrehl_Mean .* 1e6, 3)), '\n']);

fig_SlopeError = figure;
set(gca, 'FontWeight', 'bold');
hold on;
histogram(RSE_Slopes_AR1 .* 1e6, 0:0.5:20);
histogram(RSE_Slopes_CNN .* 1e6, 0:0.5:20);
histogram(RSE_Slopes_CNN_FineTuneWLS .* 1e6, 0:0.5:20);
histogram(RSE_Slopes_CNN_FineTuneStrehl .* 1e6, 0:0.5:20);
xline(pix2rad .* 1e6, 'k-', 'LineWidth', 2);
xlabel('Slope Error (\murad)');
ylabel('Counts');
legend('AR(1)', 'CNN Slope MSE', 'CNN Fine Tune WLS', ...
    'CNN Fine Tune Strehl', '1 Event Pixel', 'location', 'best');
hold off;

% exportgraphics(fig_SlopeError, 'Histogram_SlopeError_dataset04.pdf');



fig_SlopeError2 = figure;
set(gca, 'FontWeight', 'bold', 'FontSize', 12);
hold on;
histogram(RSE_Slopes_AR1 .* 1e6, 0:0.5:20);
histogram(RSE_Slopes_CNN .* 1e6, 0:0.5:20);
xline(pix2rad .* 1e6, 'k-', 'LineWidth', 2);
xlabel('Slope Error (\murad)');
ylabel('Counts');
xticks(0:2:20);
grid on;
legend('AR(1)', 'EBWFNet', '1 Event Pixel', 'location', 'best');
hold off;

% exportgraphics(fig_SlopeError2, 'Histogram_SlopeError_Dataset04_TimeAvgRem_FullApTiltRem.pdf');


%% compare OPD

OPD_FLIR = d_SlopeMSE.OPD_FLIR;
OPD_Baseline = d_SlopeMSE.OPD_baseline;
OPD_CNN_SlopeMSE = d_SlopeMSE.OPD_CNN;
OPD_CNN_FineTuneWLS = d_FineTuneWLS.OPD_CNN;
OPD_CNN_FineTuneStrehl = d_FineTuneStrehl.OPD_CNN;

% SE_OPD_Baseline = (OPD_FLIR(:) - OPD_Baseline(:)).^2;
% SE_OPD_CNN_SlopeMSE = (OPD_FLIR(:) - OPD_CNN_SlopeMSE(:)).^2;
% SE_OPD_CNN_FineTuneWLS = (OPD_FLIR(:) - OPD_CNN_FineTuneWLS(:)).^2;
% SE_OPD_CNN_FineTuneStrehl = (OPD_FLIR(:) - OPD_CNN_FineTuneStrehl(:)).^2;
% 
% % SE_OPD_Baseline_Median = median(SE_OPD_Baseline(:));
% % SE_OPD_CNN_SlopeMSE_Median = median(SE_OPD_CNN_SlopeMSE(:));
% % SE_OPD_CNN_FineTuneWLS_Median = median(SE_OPD_CNN_FineTuneWLS(:));
% % SE_OPD_CNN_FineTuneStrehl_Median = median(SE_OPD_CNN_FineTuneStrehl(:));
% 
% SE_OPD_Baseline_Mean = mean(SE_OPD_Baseline(:));
% SE_OPD_CNN_SlopeMSE_Mean = mean(SE_OPD_CNN_SlopeMSE(:));
% SE_OPD_CNN_FineTuneWLS_Mean = mean(SE_OPD_CNN_FineTuneWLS(:)); % expect this to be the smallest (best)
% SE_OPD_CNN_FineTuneStrehl_Mean = mean(SE_OPD_CNN_FineTuneStrehl(:));
% 
% fig_PhaseError = figure;
% set(gca, 'FontWeight', 'bold');
% hold on;
% histogram(SE_OPD_Baseline, 0:0.1:4);
% histogram(SE_OPD_CNN_SlopeMSE, 0:0.1:4);
% histogram(SE_OPD_CNN_FineTuneWLS, 0:0.1:4);
% histogram(SE_OPD_CNN_FineTuneStrehl, 0:0.1:4);
% xlabel('Phase Square Error (rad^{2})');
% ylabel('Counts');
% legend('Baseline', 'CNN Slope MSE', 'CNN Fine Tune WLS', ...
%     'CNN Fine Tune Strehl', 'Location', 'best');
% hold off;
% 
% % exportgraphics(fig_PhaseError, 'Histogram_PhaseError_dataset04.pdf');

%% compare Strehl

Strehl_AR1 = d_SlopeMSE.Strehl_baseline; % same in all files
Strehl_CNN_SlopeMSE = d_SlopeMSE.Strehl_CNN;
Strehl_CNN_FineTuneWLS = d_FineTuneWLS.Strehl_CNN;
Strehl_CNN_FineTuneStrehl = d_FineTuneStrehl.Strehl_CNN;

Strehl_AR1_Median = median(Strehl_AR1(:));
Strehl_CNN_SlopeMSE_Median = median(Strehl_CNN_SlopeMSE(:));
Strehl_CNN_FineTuneWLS_Median = median(Strehl_CNN_FineTuneWLS(:));
Strehl_CNN_FineTuneStrehl_Median = median(Strehl_CNN_FineTuneStrehl(:)); % expect this to be the largest (best)

Strehl_AR1_Mean = mean(Strehl_AR1(:));
Strehl_CNN_SlopeMSE_Mean = mean(Strehl_CNN_SlopeMSE(:));
Strehl_CNN_FineTuneWLS_Mean = mean(Strehl_CNN_FineTuneWLS(:));
Strehl_CNN_FineTuneStrehl_Mean = mean(Strehl_CNN_FineTuneStrehl(:)); % expect this to be the largest (best)

fprintf(1, ['Strehl Mean, AR1: ', num2str(round(Strehl_AR1_Mean, 4)), '\n']);
fprintf(1, ['Strehl Mean, CNN: ', num2str(round(Strehl_CNN_SlopeMSE_Mean, 4)), '\n']);
fprintf(1, ['Strehl Mean, WLS: ', num2str(round(Strehl_CNN_FineTuneWLS_Mean, 4)), '\n']);
fprintf(1, ['Strehl Mean, Strehl: ', num2str(round(Strehl_CNN_FineTuneStrehl_Mean, 4)), '\n']);

fig_Strehl = figure;
set(gca, 'FontWeight', 'bold');
hold on;
histogram(Strehl_AR1, 0:0.025:1);
histogram(Strehl_CNN_SlopeMSE, 0:0.025:1);
histogram(Strehl_CNN_FineTuneWLS, 0:0.025:1);
histogram(Strehl_CNN_FineTuneStrehl, 0:0.025:1);
xlabel('Approximate Strehl Ratio');
ylabel('Counts');
legend('AR(1)', 'CNN Slope MSE', 'CNN Fine Tune WLS', 'CNN Fine Tune Strehl', 'location', 'best');
hold off;

% exportgraphics(fig_Strehl, 'Histogram_Strehl_dataset04.pdf');


fig_Strehl2 = figure;
set(gca, 'FontWeight', 'bold', 'FontSize', 12);
hold on;
histogram(Strehl_AR1, 0:0.025:1);
% histogram(Strehl_CNN_SlopeMSE, 0:0.025:1);
histogram(Strehl_CNN_FineTuneStrehl, 0:0.025:1);
xlabel('Approximate Strehl Ratio');
ylabel('Counts');
xticks(0:0.1:1);
grid on;
legend('AR(1)', 'EBWFNet', 'location', 'northwest');
hold off;

% exportgraphics(fig_Strehl2, 'Histogram_Strehl_Dataset04_TimeAvgRem_FullApTiltRem.pdf');



[f_Strehl_Baseline, x_Strehl_Baseline] = ecdf(Strehl_AR1(:));
[f_Strehl_CNN_SlopeMSE, x_Strehl_CNN_SlopeMSE] = ecdf(Strehl_CNN_SlopeMSE(:));
[f_Strehl_CNN_FineTuneWLS, x_Strehl_CNN_FineTuneWLS] = ecdf(Strehl_CNN_FineTuneWLS(:));
[f_Strehl_CNN_FineTuneStrehl, x_Strehl_CNN_FineTuneStrehl] = ecdf(Strehl_CNN_FineTuneStrehl(:));

% figure;
% set(gca, 'FontWeight', 'bold', 'FontSize', 12);
% hold on;
% plot(x_Strehl_Baseline, f_Strehl_Baseline, '-', 'LineWidth', 2);
% plot(x_Strehl_CNN_SlopeMSE, f_Strehl_CNN_SlopeMSE, '-', 'LineWidth', 2);
% plot(x_Strehl_CNN_FineTuneWLS, f_Strehl_CNN_FineTuneWLS, '-', 'LineWidth', 2);
% plot(x_Strehl_CNN_FineTuneStrehl, f_Strehl_CNN_FineTuneStrehl, '-', 'LineWidth', 2);
% xlim([0, 1]);
% ylim([0, 1]);
% xlabel('Approximate Strehl Ratio');
% ylabel('CDF, P[X \leq x]');
% grid on;
% legend('Avg. Position', 'CNN Slope MSE', 'CNN Fine Tune WLS', 'CNN Fine Tune Strehl', 'location', 'northwest');
% hold off;


%%

Recon = d_SlopeMSE.Recon; % reconstructor

[U,S,V] = svd(Recon); % SVD of reconstructor

[numOPDs, numSlopes] = size(S);

error_slopeSpace_Baseline = V' * (slopes_Frame - slopes_AR1);
error_slopeSpace_CNN_SlopeMSE = V' * (slopes_Frame - slopes_CNN);
error_slopeSpace_CNN_FineTuneWLS = V' * (slopes_Frame - slopes_CNN_FineTuneWLS);
error_slopeSpace_CNN_FineTuneStrehl = V' * (slopes_Frame - slopes_CNN_FineTuneStrehl);

error_reconSpace_Baseline = S * error_slopeSpace_Baseline;
error_reconSpace_CNN_SlopeMSE = S * error_slopeSpace_CNN_SlopeMSE;
error_reconSpace_CNN_FineTuneWLS = S * error_slopeSpace_CNN_FineTuneWLS;
error_reconSpace_CNN_FineTuneStrehl = S * error_slopeSpace_CNN_FineTuneStrehl;

var_error_slopeSpace_Baseline = var(error_slopeSpace_Baseline, [], 2);
var_error_slopeSpace_CNN_SlopeMSE = var(error_slopeSpace_CNN_SlopeMSE, [], 2);
var_error_slopeSpace_CNN_FineTuneWLS = var(error_slopeSpace_CNN_FineTuneWLS, [], 2);
var_error_slopeSpace_CNN_FineTuneStrehl = var(error_slopeSpace_CNN_FineTuneStrehl, [], 2);

var_error_reconSpace_Baseline = var(error_reconSpace_Baseline, [], 2);
var_error_reconSpace_CNN_SlopeMSE = var(error_reconSpace_CNN_SlopeMSE, [], 2);
var_error_reconSpace_CNN_FineTuneWLS = var(error_reconSpace_CNN_FineTuneWLS, [], 2);
var_error_reconSpace_CNN_FineTuneStrehl = var(error_reconSpace_CNN_FineTuneStrehl, [], 2);

sum_var_error_slopeSpace_Baseline = sum(var_error_slopeSpace_Baseline);
sum_var_error_slopeSpace_CNN_SlopeMSE = sum(var_error_slopeSpace_CNN_SlopeMSE); % expect this to be smallest?
sum_var_error_slopeSpace_CNN_FineTuneWLS = sum(var_error_slopeSpace_CNN_FineTuneWLS); % expect this to be smallest??
sum_var_error_slopeSpace_CNN_FineTuneStrehl = sum(var_error_slopeSpace_CNN_FineTuneStrehl);

sum_var_error_reconSpace_Baseline = sum(var_error_reconSpace_Baseline);
sum_var_error_reconSpace_CNN_SlopeMSE = sum(var_error_reconSpace_CNN_SlopeMSE);
sum_var_error_reconSpace_CNN_FineTuneWLS = sum(var_error_reconSpace_CNN_FineTuneWLS);
sum_var_error_reconSpace_CNN_FineTuneStrehl = sum(var_error_reconSpace_CNN_FineTuneStrehl); % expect this to be smallest

% sum_var_error_slopeSpace_Baseline = sum(var(error_slopeSpace_Baseline(:)));
% sum_var_error_slopeSpace_CNN_SlopeMSE = sum(var(error_slopeSpace_CNN_SlopeMSE(:)));
% sum_var_error_slopeSpace_CNN_FineTuneWLS = sum(var(error_slopeSpace_CNN_FineTuneWLS(:))); % expect this to be smallest
% sum_var_error_slopeSpace_CNN_FineTuneStrehl = sum(var(error_slopeSpace_CNN_FineTuneStrehl(:)));
% 
% sum_var_error_reconSpace_Baseline = sum(var(error_reconSpace_Baseline(:)));
% sum_var_error_reconSpace_CNN_SlopeMSE = sum(var(error_reconSpace_CNN_SlopeMSE(:)));
% sum_var_error_reconSpace_CNN_FineTuneWLS = sum(var(error_reconSpace_CNN_FineTuneWLS(:)));
% sum_var_error_reconSpace_CNN_FineTuneStrehl = sum(var(error_reconSpace_CNN_FineTuneStrehl(:))); % expect this to be smallest

fig1 = figure;
set(gca, 'FontWeight', 'bold');
hold on;
plot(1:numOPDs, diag(S), 'LineWidth', 2);
xlim([0, numOPDs+1]);
xlabel('Singular Vector');
ylabel('Singular Value');
title('Reconstructor Singular Values');
hold off;

% exportgraphics(fig1, 'ReconstructorSingularValues.pdf');

fig2 = figure;
set(gca, 'YScale', 'log', 'FontWeight', 'bold');
hold on;
plot(1:numSlopes, var_error_slopeSpace_Baseline, 'LineWidth', 2);
plot(1:numSlopes, var_error_slopeSpace_CNN_SlopeMSE, 'LineWidth', 2);
plot(1:numSlopes, var_error_slopeSpace_CNN_FineTuneWLS, 'LineWidth', 2);
plot(1:numSlopes, var_error_slopeSpace_CNN_FineTuneStrehl, 'LineWidth', 2);
xlim([0, numSlopes+1]);
% xlim([0, 15]);
xlabel('Singular Vector');
ylabel('Error Variance');
title('Slope Space');
% legend('Baseline', 'CNN Slope MSE', 'CNN Fine Tune WLS', 'CNN Fine Tune Strehl', 'location', 'northwest');
legend('Baseline', 'CNN Slope MSE', 'CNN Fine Tune WLS', 'CNN Fine Tune Strehl', 'location', 'northeast');
hold off;

% % exportgraphics(fig2, 'ErrorVariance_SlopeSpace_Zoom.pdf');
% exportgraphics(fig2, 'ErrorVariance_SlopeSpace.pdf');

fig3 = figure;
set(gca, 'YScale', 'log', 'FontWeight', 'bold');
hold on;
plot(1:numOPDs, var_error_reconSpace_Baseline, 'LineWidth', 2);
plot(1:numOPDs, var_error_reconSpace_CNN_SlopeMSE, 'LineWidth', 2);
plot(1:numOPDs, var_error_reconSpace_CNN_FineTuneWLS, 'LineWidth', 2);
plot(1:numOPDs, var_error_reconSpace_CNN_FineTuneStrehl, 'LineWidth', 2);
xlim([0, numOPDs+1]);
% xlim([0, 15]);
ylim([1e-18, 1e-12]);
xlabel('Singular Vector');
ylabel('Error Variance');
title('Reconstruction Space');
% legend('Baseline', 'CNN Slope MSE', 'CNN Fine Tune WLS', 'CNN Fine Tune Strehl', 'location', 'northeast');
legend('Baseline', 'CNN Slope MSE', 'CNN Fine Tune WLS', 'CNN Fine Tune Strehl', 'location', 'southwest');
hold off;

% % exportgraphics(fig3, 'ErrorVariance_ReconstructionSpace_Zoom.pdf');
% exportgraphics(fig3, 'ErrorVariance_ReconstructionSpace.pdf');

%% compute per-frame structure function

numFrames = size(OPD_FLIR, 2);

xact_custom = d_SlopeMSE.xact_custom;
yact_custom = d_SlopeMSE.yact_custom;
dx_act = d_SlopeMSE.dx_act;

[opdrec2d_FLIR_tmp, x1, y1] = actmap(OPD_FLIR(:,1), xact_custom, yact_custom, [], [], 20);

mask = ~isnan(opdrec2d_FLIR_tmp);

for ii=1:numFrames
    [opdrec2d_FLIR_tmp, ~, ~] = actmap(OPD_FLIR(:,ii), xact_custom, yact_custom, [], [], 20);
    [opdrec2d_Baseline_tmp, ~, ~] = actmap(OPD_Baseline(:,ii), xact_custom, yact_custom, [], [], 20);
    [opdrec2d_CNN_SlopeMSE_tmp, ~, ~] = actmap(OPD_CNN_SlopeMSE(:,ii), xact_custom, yact_custom, [], [], 20);
    [opdrec2d_CNN_FineTuneWLS_tmp, ~, ~] = actmap(OPD_CNN_FineTuneWLS(:,ii), xact_custom, yact_custom, [], [], 20);
    [opdrec2d_CNN_FineTuneStrehl_tmp, ~, ~] = actmap(OPD_CNN_FineTuneStrehl(:,ii), xact_custom, yact_custom, [], [], 20);

    opdrec2d_FLIR(:,:,ii) = opdrec2d_FLIR_tmp;
    opdrec2d_Baseline(:,:,ii) = opdrec2d_Baseline_tmp;
    opdrec2d_CNN_SlopeMSE(:,:,ii) = opdrec2d_CNN_SlopeMSE_tmp;
    opdrec2d_CNN_FineTuneWLS(:,:,ii) = opdrec2d_CNN_FineTuneWLS_tmp;
    opdrec2d_CNN_FineTuneStrehl(:,:,ii) = opdrec2d_CNN_FineTuneStrehl_tmp;

    opdrec2d_FLIR_tmp(~mask) = 0;
    opdrec2d_Baseline_tmp(~mask) = 0;
    opdrec2d_CNN_SlopeMSE_tmp(~mask) = 0;
    opdrec2d_CNN_FineTuneWLS_tmp(~mask) = 0;
    opdrec2d_CNN_FineTuneStrehl_tmp(~mask) = 0;

    [D_FLIR_tmp, ~] = str_fcn2_ft_new(opdrec2d_FLIR_tmp, mask, dx_act);
    [D_Baseline_tmp, ~] = str_fcn2_ft_new(opdrec2d_Baseline_tmp, mask, dx_act);
    [D_CNN_SlopeMSE_tmp, ~] = str_fcn2_ft_new(opdrec2d_CNN_SlopeMSE_tmp, mask, dx_act);
    [D_CNN_FineTuneWLS_tmp, ~] = str_fcn2_ft_new(opdrec2d_CNN_FineTuneWLS_tmp, mask, dx_act);
    [D_CNN_FineTuneStrehl_tmp, ~] = str_fcn2_ft_new(opdrec2d_CNN_FineTuneStrehl_tmp, mask, dx_act);

    D_FLIR(:,:,ii) = D_FLIR_tmp;
    D_Baseline(:,:,ii) = D_Baseline_tmp;
    D_CNN_SlopeMSE(:,:,ii) = D_CNN_SlopeMSE_tmp;
    D_CNN_FineTuneWLS(:,:,ii) = D_CNN_FineTuneWLS_tmp;
    D_CNN_FineTuneStrehl(:,:,ii) = D_CNN_FineTuneStrehl_tmp;
end

clear('D_FLIR_tmp', 'D_Baseline_tmp', 'D_CNN_SlopeMSE_tmp', 'D_CNN_FineTuneWLS_tmp', ...
    'D_CNN_FineTuneStrehl_tmp', 'opdrec2d_FLIR_tmp', ...
    'opdrec2d_Baseline_tmp', 'opdrec2d_CNN_SlopeMSE_tmp', ...
    'opdrec2d_CNN_FineTuneWLS_tmp', 'opdrec2d_CNN_FineTuneStrehl_tmp');

cmin = min([D_FLIR(:,:,1); D_Baseline(:,:,1); D_CNN_SlopeMSE(:,:,1)], [], 'all');
cmax = max([D_FLIR(:,:,1); D_Baseline(:,:,1); D_CNN_SlopeMSE(:,:,1)], [], 'all');

% figure('Position', [100, 100, 1200, 500]);
% subplot(1, 3, 1);
% imagesc(x1, y1, D_FLIR(:,:,1));
% axis image xy square;
% colorbar; clim([cmin, cmax]);
% title('FLIR');
% 
% subplot(1, 3, 2);
% imagesc(x1, y1, D_Baseline(:,:,1));
% axis image xy square;
% colorbar; clim([cmin, cmax]);
% title('Baseline');
% 
% subplot(1, 3, 3);
% imagesc(x1, y1, D_CNN_SlopeMSE(:,:,1));
% axis image xy square;
% colorbar; clim([cmin, cmax]);
% title('CNN Slope MSE');

%% analyze 2D structure functions

% time-average structure function
D_FLIR_Mean = mean(D_FLIR, 3);
D_Baseline_Mean = mean(D_Baseline, 3);
D_CNN_SlopeMSE_Mean = mean(D_CNN_SlopeMSE, 3);
D_CNN_FineTuneWLS_Mean = mean(D_CNN_FineTuneWLS, 3);
D_CNN_FineTuneStrehl_Mean = mean(D_CNN_FineTuneStrehl, 3);

% per-frame squared error
D_Baseline_SE = (D_FLIR - D_Baseline).^2;
D_CNN_SlopeMSE_SE = (D_FLIR - D_CNN_SlopeMSE).^2;
D_CNN_FineTuneWLS_SE = (D_FLIR - D_CNN_FineTuneWLS).^2;
D_CNN_FineTuneStrehl_SE = (D_FLIR - D_CNN_FineTuneStrehl).^2;

% time-average of per-frame squared error
D_Baseline_SE_Mean = mean(D_Baseline_SE, 3);
D_CNN_SlopeMSE_SE_Mean = mean(D_CNN_SlopeMSE_SE, 3);
D_CNN_FineTuneWLS_SE_Mean = mean(D_CNN_FineTuneWLS_SE, 3);
D_CNN_FineTuneStrehl_SE_Mean = mean(D_CNN_FineTuneStrehl_SE, 3);

% time-variance of per-frame squared error
D_Baseline_SE_Var = var(D_Baseline_SE, [], 3);
D_CNN_SlopeMSE_SE_Var = var(D_CNN_SlopeMSE_SE, [], 3);
D_CNN_FineTuneWLS_SE_Var = var(D_CNN_FineTuneWLS_SE, [], 3);
D_CNN_FineTuneStrehl_SE_Var = var(D_CNN_FineTuneStrehl_SE, [], 3);

cmin = 0;
% cmin = min([D_Baseline_SE_Mean(:); D_CNN_SlopeMSE_SE_Mean(:); D_CNN_FineTuneWLS_SE_Mean(:); D_CNN_FineTuneStrehl_SE_Mean(:)]);
cmax = max([D_Baseline_SE_Mean(:); D_CNN_SlopeMSE_SE_Mean(:); D_CNN_FineTuneWLS_SE_Mean(:); D_CNN_FineTuneStrehl_SE_Mean(:)]);

fig_D_Baseline_SE_Mean = figure;
% imagesc(x1, y1, D_Baseline_SE_Mean);
imagesc(y1, x1, D_Baseline_SE_Mean');
axis image xy square;
colormap(jet);
colorbar; clim([cmin, cmax]);
xlim([-0.2, 0.2]);
ylim([-0.2, 0.2]);
xlabel('x (m)');
ylabel('y (m)');
% title('Mean-Square-Error: Average Position Tracker');
set(gca, 'FontWeight', 'bold');

% exportgraphics(fig_D_Baseline_SE_Mean, 'Figure_StructureFunction_AvgPosTracker_TimeAvgRmv_FullApTiltInc_SE_Mean.pdf');


fig_D_CNN_SlopeMSE_SE_Mean = figure;
% imagesc(x1, y1, D_CNN_SlopeMSE_SE_Mean);
imagesc(y1, x1, D_CNN_SlopeMSE_SE_Mean');
axis image xy square;
colormap(jet);
colorbar; clim([cmin, cmax]);
xlim([-0.2, 0.2]);
ylim([-0.2, 0.2]);
xlabel('x (m)');
ylabel('y (m)');
% title('Mean-Square-Error: CNN Slope MSE');
set(gca, 'FontWeight', 'bold');

% exportgraphics(fig_D_CNN_SlopeMSE_SE_Mean, 'Figure_StructureFunction_CNN_SlopeMSE_TimeAvgRmv_FullApTiltInc_SE_Mean.pdf');


% fig_D_CNN_FineTuneWLS_SE_Mean = figure;
% % imagesc(x1, y1, D_CNN_FineTuneWLS_SE_Mean);
% imagesc(y1, x1, D_CNN_FineTuneWLS_SE_Mean');
% axis image xy square;
% colormap(jet);
% colorbar; clim([cmin, cmax]);
% xlim([-0.2, 0.2]);
% ylim([-0.2, 0.2]);
% xlabel('x (m)');
% ylabel('y (m)');
% title('Mean-Square-Error: CNN Fine Tune WLS');
% set(gca, 'FontWeight', 'bold');
% 
% % exportgraphics(fig_D_CNN_FineTuneWLS_SE_Mean, 'Figure_StructureFunction_CNN_FineTuneWLS_SE_Mean.pdf');


% fig_D_CNN_FineTuneStrehl_SE_Mean = figure;
% % imagesc(x1, y1, D_CNN_FineTuneStrehl_SE_Mean);
% imagesc(y1, x1, D_CNN_FineTuneStrehl_SE_Mean');
% axis image xy square;
% colormap(jet);
% colorbar; clim([cmin, cmax]);
% xlim([-0.2, 0.2]);
% ylim([-0.2, 0.2]);
% xlabel('x (m)');
% ylabel('y (m)');
% title('Mean-Square-Error: CNN Fine Tune Strehl');
% set(gca, 'FontWeight', 'bold');
% 
% % exportgraphics(fig_D_CNN_FineTuneStrehl_SE_Mean, 'Figure_StructureFunction_CNN_FineTuneStrehl_SE_Mean.pdf');


cmin = 0;
% cmin = min([D_Baseline_SE_Var(:); D_CNN_SlopeMSE_SE_Var(:); D_CNN_FineTuneWLS_SE_Var(:); D_CNN_FineTuneStrehl_SE_Var(:)]);
cmax = max([D_Baseline_SE_Var(:); D_CNN_SlopeMSE_SE_Var(:); D_CNN_FineTuneWLS_SE_Var(:); D_CNN_FineTuneStrehl_SE_Var(:)]);

fig_D_Baseline_SE_Var = figure;
% imagesc(x1, y1, D_Baseline_SE_Var);
imagesc(y1, x1, D_Baseline_SE_Var');
axis image xy square;
colormap(jet);
colorbar; clim([cmin, cmax]);
xlim([-0.2, 0.2]);
ylim([-0.2, 0.2]);
xlabel('x (m)');
ylabel('y (m)');
title('Variance-Square-Error: Average Position Tracker');
set(gca, 'FontWeight', 'bold');

% exportgraphics(fig_D_Baseline_SE_Var, 'Figure_StructureFunction_Baseline_SE_Var.pdf');


fig_D_CNN_SlopeMSE_SE_Var = figure;
% imagesc(x1, y1, D_CNN_SlopeMSE_SE_Var);
imagesc(y1, x1, D_CNN_SlopeMSE_SE_Var');
axis image xy square;
colormap(jet);
colorbar; clim([cmin, cmax]);
xlim([-0.2, 0.2]);
ylim([-0.2, 0.2]);
xlabel('x (m)');
ylabel('y (m)');
title('Variance-Square-Error: CNN Slope MSE');
set(gca, 'FontWeight', 'bold');

% exportgraphics(fig_D_CNN_SlopeMSE_SE_Var, 'Figure_StructureFunction_CNN_SlopeMSE_SE_Var.pdf');


% fig_D_CNN_FineTuneWLS_SE_Var = figure;
% % imagesc(x1, y1, D_CNN_FineTuneWLS_SE_Var);
% imagesc(y1, x1, D_CNN_FineTuneWLS_SE_Var');
% axis image xy square;
% colormap(jet);
% colorbar; clim([cmin, cmax]);
% xlim([-0.2, 0.2]);
% ylim([-0.2, 0.2]);
% xlabel('x (m)');
% ylabel('y (m)');
% title('Variance-Square-Error: CNN Fine Tune WLS');
% set(gca, 'FontWeight', 'bold');
% 
% % exportgraphics(fig_D_CNN_FineTuneWLS_SE_Var, 'Figure_StructureFunction_CNN_FineTuneWLS_SE_Var.pdf');


% fig_D_CNN_FineTuneStrehl_SE_Var = figure;
% % imagesc(x1, y1, D_CNN_FineTuneStrehl_SE_Var);
% imagesc(y1, x1, D_CNN_FineTuneStrehl_SE_Var');
% axis image xy square;
% colormap(jet);
% colorbar; clim([cmin, cmax]);
% xlim([-0.2, 0.2]);
% ylim([-0.2, 0.2]);
% xlabel('x (m)');
% ylabel('y (m)');
% title('Variance-Square-Error: CNN Fine Tune Strehl');
% set(gca, 'FontWeight', 'bold');
% 
% % exportgraphics(fig_D_CNN_FineTuneStrehl_SE_Var, 'Figure_StructureFunction_CNN_FineTuneStrehl_SE_Var.pdf');

%% azimuthal-average of structure function

[r_FLIR, dataAvg_FLIR] = azimuthalAvg(x1, y1, D_FLIR_Mean');
[r_Baseline, dataAvg_Baseline] = azimuthalAvg(x1, y1, D_Baseline_Mean');
[r_CNN_SlopeMSE, dataAvg_CNN_SlopeMSE] = azimuthalAvg(x1, y1, D_CNN_SlopeMSE_Mean');
[r_CNN_FineTuneWLS, dataAvg_CNN_FineTuneWLS] = azimuthalAvg(x1, y1, D_CNN_FineTuneWLS_Mean');
[r_CNN_FineTuneStrehl, dataAvg_CNN_FineTuneStrehl] = azimuthalAvg(x1, y1, D_CNN_FineTuneStrehl_Mean');

idx_use = find(r_FLIR >= 0.02 & r_FLIR <= 0.12);
r0_theory_tmp = 0.01:0.001:0.1; % 1cm to 10cm
for ii = 1:length(r0_theory_tmp)
    D_theory_tmp = 6.88 .* (r_FLIR ./ r0_theory_tmp(ii)).^(5/3);

    D_FLIR_Theory_MSE(ii) = mean((D_theory_tmp(idx_use) - dataAvg_FLIR(idx_use)).^2);
    D_Baseline_Theory_MSE(ii) = mean((D_theory_tmp(idx_use) - dataAvg_Baseline(idx_use)).^2);
    D_CNN_SlopeMSE_Theory_MSE(ii) = mean((D_theory_tmp(idx_use) - dataAvg_CNN_SlopeMSE(idx_use)).^2);
    D_CNN_FineTuneWLS_Theory_MSE(ii) = mean((D_theory_tmp(idx_use) - dataAvg_CNN_FineTuneWLS(idx_use)).^2);
    D_CNN_FineTuneStrehl_Theory_MSE(ii) = mean((D_theory_tmp(idx_use) - dataAvg_CNN_FineTuneStrehl(idx_use)).^2);
end

[~, I_min] = min(D_FLIR_Theory_MSE);
r0_Theory = r0_theory_tmp(I_min);
D_Theory = 6.88 .* (r_FLIR ./ r0_Theory).^(5/3);

% figure;
% plot(r0_theory_tmp.*100, D_FLIR_Theory_MSE, 'LineWidth', 2);
% hold on;
% plot(r0_theory_tmp.*100, D_Baseline_Theory_MSE, 'LineWidth', 2);
% plot(r0_theory_tmp.*100, D_CNN_SlopeMSE_Theory_MSE, 'LineWidth', 2);
% plot(r0_theory_tmp.*100, D_CNN_FineTuneWLS_Theory_MSE, 'LineWidth', 2);
% plot(r0_theory_tmp.*100, D_CNN_FineTuneStrehl_Theory_MSE, 'LineWidth', 2);
% xlabel('r_{0} (cm)');
% ylabel('MSE');
% hold off;

D_FLIR_Baseline_MSE = mean((dataAvg_FLIR(idx_use) - dataAvg_Baseline(idx_use)).^2);
D_FLIR_CNN_SlopeMSE_MSE = mean((dataAvg_FLIR(idx_use) - dataAvg_CNN_SlopeMSE(idx_use)).^2);
D_FLIR_CNN_FineTuneWLS_MSE = mean((dataAvg_FLIR(idx_use) - dataAvg_CNN_FineTuneWLS(idx_use)).^2);
D_FLIR_CNN_FineTuneStrehl_MSE = mean((dataAvg_FLIR(idx_use) - dataAvg_CNN_FineTuneStrehl(idx_use)).^2);


% fig_StructureFunction_AzAvg = figure;
% plot(r_FLIR.*100, dataAvg_FLIR, '-', 'LineWidth', 2);
% hold on;
% plot(r_Baseline.*100, dataAvg_Baseline, '-', 'LineWidth', 2);
% plot(r_CNN_SlopeMSE.*100, dataAvg_CNN_SlopeMSE, '-', 'LineWidth', 2);
% plot(r_CNN_FineTuneWLS.*100, dataAvg_CNN_FineTuneWLS, '-', 'LineWidth', 2);
% plot(r_CNN_FineTuneStrehl.*100, dataAvg_CNN_FineTuneStrehl, '-', 'LineWidth', 2);
% plot(r_FLIR.*100, D_Theory, 'k-', 'LineWidth', 2);
% xlim([0.02, 0.15].*100);
% set(gca, 'XScale', 'log', 'YScale', 'log', 'FontWeight', 'bold');
% xlabel('Radial Separation, r (cm)');
% ylabel('Phase Structure Function D_{\phi}(r)');
% legend('FLIR', 'Avg. Position', 'CNN Slope MSE', 'CNN Fine Tune WLS', 'CNN Fine Tune Strehl', ...
%     ['r_{0} = ', num2str(r0_Theory * 100), ' cm'], 'location', 'northwest');
% hold off;
% 
% % exportgraphics(fig_StructureFunction_AzAvg, 'Figure_StructureFunction_AzAvg.pdf');


fig_StructureFunction_AzAvg = figure;
plot(r_FLIR.*100, dataAvg_FLIR, 'k-', 'LineWidth', 2);
hold on;
plot(r_Baseline.*100, dataAvg_Baseline, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
plot(r_CNN_SlopeMSE.*100, dataAvg_CNN_SlopeMSE, '--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
plot(r_FLIR.*100, D_Theory, '-', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2);
xlim([0.02, 0.15].*100);
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontWeight', 'bold', 'FontSize', 12);
xlabel('Radial Separation, r (cm)');
ylabel('Phase Structure Function D_{\phi}(r)');
grid on;
legend('Frame', 'AR(1)', 'EBWFNet', ...
    ['r_{0} = ', num2str(r0_Theory * 100), ' cm'], 'location', 'northwest');
hold off;

% exportgraphics(fig_StructureFunction_AzAvg, 'Figure_StructureFunction_AzAvg_TimeAvgRmv_FullApTiltInc.pdf');
