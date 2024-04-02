clear; close all; clc;
% 
% % x = 1:1000;
% x = linspace(-pi, pi, 1000);
% y = sin(x);
% 
% figure;
% % plot(x, y);
% plot(1:1000, y);

%% generate autoregressive process
% https://www.mathworks.com/help/signal/ug/linear-prediction-and-autoregressive-modeling.html

b = fir1(1024, .0125, 'low');
[d, p0] = lpc(b, 7);

rng(0, 'twister'); % Allow reproduction of exact experiment
u = sqrt(p0)*randn(1000,1); % White Gaussian noise with variance p0

x = filter(1, d, u);

% only consider first 500 samples
x2 = x(1:500);

%%

% normalize intensity
I = 8 * (x2 - min(x2)) / (max(x2) - min(x2)) + 1.65; % +0.99 because weird stuff happens with +1

I_floor = floor(I);

%%

idx_pos = [];
idx_neg = [];

I_level = round(I(1));
for ii = 2:length(I)
    I_now = I(ii);
    if I_now >= I_level + 1
        idx_pos = [idx_pos, ii];
        I_level = I_level + 1;
    elseif I_now <= I_level - 1
        idx_neg = [idx_neg, ii];
        I_level = I_level - 1;
    end
end

fig_Intensity = figure('Position', [600, 600, 1000, 200]);
for ii = 1:10
    plot([1, 500], [ii, ii], '-', 'LineWidth', 2, 'Color', [1, 1, 1] .* 0.8);
    hold on;
end
h1 = plot(I, 'k-', 'LineWidth', 2);
h2 = stem(idx_pos, I_floor(idx_pos), 'Marker', 'None', 'LineWidth', 2, 'Color', 'blue');
h3 = stem(idx_neg, I_floor(idx_neg)+1, 'Marker', 'None', 'LineWidth', 2, 'Color', 'red');
yticks(0:10);
ylim([0, 10]);
xlim([0, 450]);
xlabel('Time (ms)');
ylabel('Levels (Log Scale)');
grid on;
legend([h1, h2, h3], {'Intensity', 'Positive Event', 'Negative Event'}, 'location', 'North');
set(gca, 'FontWeight', 'bold', 'FontSize', 12);
hold off;

fig_Events = figure('Position', [600, 200, 1000, 200]);
plot([1, 500], [0, 0], '-', 'LineWidth', 2, 'Color', [1, 1, 1] .* 0.8);
hold on;
h1 = stem(idx_pos, ones(1, length(idx_pos)), 'filled', 'Color', 'blue', 'LineWidth', 2);
h2 = stem(idx_neg, -1 .* ones(1, length(idx_neg)), 'filled', 'Color', 'red', 'LineWidth', 2);
ylim([-1.25, 1.25]);
xlim([0, 450]);
yticks(-1:1);
xlabel('Time (ms)');
ylabel('Polarity');
grid on;
legend([h1, h2], {'Positive Event', 'Negative Event'});
set(gca, 'FontWeight', 'bold', 'FontSize', 12);
hold off;


% exportgraphics(fig_Intensity, 'Figure_Synthetic_Intensity_Events.pdf');
% exportgraphics(fig_Events, 'Figure_Synthetic_Events_Stem.pdf');

%% another (wrong?) implementation
% x_plot = [1,500];
% y_plot = 1:10;
% 
% figure('Position', [600, 600, 1000, 300]);
% % yline(1:10, '-', 'LineWidth', 2, 'Color', [0.5, 0.5, 0.5]);
% for ii = 1:10
%     plot([1, 500], [ii, ii], '-', 'LineWidth', 2, 'Color', [1, 1, 1] .* 0.8);
%     hold on;
% end
% plot(I, 'k-', 'LineWidth', 2);
% yticks(0:11);
% ylim([0, 11]);
% xlabel('Time (ms)');
% ylabel('Levels (Log Scale)');
% set(gca, 'FontWeight', 'bold');
% hold off;
% 
% idx_pos = find(diff(I_floor) > 0) + 1;
% idx_neg = find(diff(I_floor) < 0) + 1;
% 
% figure('Position', [600, 200, 1000, 300]);
% plot([1, 500], [0, 0], '-', 'LineWidth', 2, 'Color', [1, 1, 1] .* 0.8);
% hold on;
% stem(idx_pos, ones(1, length(idx_pos)), 'filled', 'Color', 'blue', 'LineWidth', 2);
% stem(idx_neg, -1 .* ones(1, length(idx_neg)), 'filled', 'Color', 'red', 'LineWidth', 2);
% ylim([-1.25, 1.25]);
% yticks(-1:1);
% 
% figure('Position', [600, 600, 1000, 300]);
% for ii = 1:10
%     plot([1, 500], [ii, ii], '-', 'LineWidth', 2, 'Color', [1, 1, 1] .* 0.8);
%     hold on;
% end
% plot(I, 'k-', 'LineWidth', 2);
% stem(idx_pos, I_floor(idx_pos), 'Marker', 'None', 'LineWidth', 2, 'BaseValue', 0);
% stem(idx_neg, I_floor(idx_neg)+1, 'Marker', 'None', 'LineWidth', 2, 'BaseValue', 0);
% yticks(0:11);
% ylim([0, 11]);
% xlabel('Time (ms)');
% ylabel('Levels (Log Scale)');
% set(gca, 'FontWeight', 'bold');
% hold off;