clear; close all; clc;
addpath('..\..\Filtering');

% Script that processes frame (FLIR) data by reading in the data, removing
% background for threshold center-of-gravity (TCoG) spot position
% estimation. Data are saved for further processing.
% Author: Mitchell Grose, University of Dayton, 2023

% % event tracking threshold
track_thresh_ev = 0.25;

evFrames_FPS = 100; % desired event frames FPS
FLIR_FPS = 100; % FLIR frames captured at this FPS

% subaperture sizes
subap_sz_ev = 54; % Prophesee VGA sensor

% hardware specs
M1 = -30e-3 / 1250e-3; % front 4-f system
M2 = -100e-3 / 40e-3; % relay lens
f_MLA = 24e-3; % focal length of MLA

% pixel pitch of each sensor
pitch_EV = 15e-6; % Prophesee VGA, 15 micrometer pitch
pitch_FLIR = 3.45e-6; % FLIR, 3.45 micrometer pitch

% FLIR data processing thresholds
thresh_reference = 0.4;
thresh_global = 0.20;

% subaperture sizes
subap_sz_flir = 200;

%% data to use

date_flag = '2023-07-30';
dataset_flag = 4;

fpath_ev = ['D:\Data\Dual_WFS\', date_flag, '_Eastwood\DualSync\Event\dataset', num2str(dataset_flag, '%02d')];
fpath_flir = ['D:\Data\Dual_WFS\', date_flag, '_Eastwood\DualSync\FLIR\dataset', num2str(dataset_flag, '%02d')];

dd_flir = dir(fullfile(fpath_flir, '*.tiff'));
% dd_flir = dd_flir(1:900);

firstframe_keep = 1;

% numFramesLoad = length(dd_flir);
numFramesLoad = 600;

% load images and format
imgs = zeros([1536, 2048, numFramesLoad], 'uint8');
% flirFrame_ref = zeros([1536, 2048]);

fprintf(1, 'reading FLIR images... ');
for ii=1:numFramesLoad
    if rem(ii, 10) == 0
        fprintf(1, ['Loading frame ', num2str(ii), ' of ', num2str(numFramesLoad), '.\n']);
    end
    % fliplr because 45 degree beamsplitter flips the image left-right
    % flirFrame_ref = flirFrame_ref + fliplr(double(imread(fullfile(dd_FLIR(ii).folder, dd_FLIR(ii).name))));
    imgs(:,:,ii) = fliplr(imread(fullfile(dd_flir(ii).folder, dd_flir(ii).name)));
end

% flirFrame_ref = flirFrame_ref / numFramesLoad;
flirFrame_ref = mean(imgs, 3);

flirFrame_ref = flirFrame_ref - max(flirFrame_ref(:)) .* thresh_reference;
flirFrame_ref(flirFrame_ref < 0) = 0;

figure; imagesc(flirFrame_ref); axis image xy; colormap(gray); colorbar; drawnow;

%% apply global threshold

% apply global threshold to images
tic; fprintf(1, 'Filtering images with global threshold... ');
imgs_thresh = single(imgs);
clear('imgs');

for ii=1:size(imgs_thresh, 3)
    if rem(ii, 10) == 0
        fprintf(1, ['Thresholding frame ', num2str(ii), ' of ', num2str(numFramesLoad), '.\n']);
    end
    img_tmp = imgs_thresh(:,:,ii);
    img_tmp = img_tmp - (max(img_tmp(:)) .* thresh_global);
    img_tmp(img_tmp<0) = 0;
    imgs_thresh(:,:,ii) = img_tmp;

    % imgs_thresh(:,:,ii) = imgs_thresh(:,:,ii) - max(imgs_thresh(:,:,ii), [], 'all').*thresh_global;
end
% imgs_thresh(imgs_thresh<0) = 0;

clear('img_tmp');
fprintf(1, 'Done. '); toc;

fig1 = figure;
hold on;
for ii = size(imgs_thresh, 3)-400:size(imgs_thresh, 3)
    imagesc(imgs_thresh(:,:,ii));
    axis image xy;
    colormap(gray);
    colorbar;
    title(['Frame ', num2str(ii)]);

    drawnow;
    clf(fig1);
end
close(fig1);

%% load reference positions

load(['./spot_ref_centers_20230730/ref_pos_flir_dataset', num2str(dataset_flag, '%02d'), '.mat']);

% define rectangles for subaperture positions
r_pos_flir = [xc_ref_flir-subap_sz_flir/2, yc_ref_flir-subap_sz_flir/2,...
    ones(size(yc_ref_flir)).*subap_sz_flir, ones(size(yc_ref_flir)).*subap_sz_flir];

flirFrame_ref_padded = padarray(flirFrame_ref, [50, 50], 0, 'both');
r_pos_flir_padded = r_pos_flir;
r_pos_flir_padded(:,1:2) = r_pos_flir_padded(:,1:2) + 50;

% plot reference image with subaperture positions
fig1 = figure;
% imagesc(flirFrame_ref);
imagesc(flirFrame_ref_padded);
hold on;
% plot(xc_ref_flir, yc_ref_flir, 'r+', 'LineWidth', 2, 'MarkerSize', 6);
for ii=1:length(xc_ref_flir)
    % rectangle('Position', r_pos_flir(ii,:), 'EdgeColor', 'g', 'LineWidth', 3);
    rectangle('Position', r_pos_flir_padded(ii,:), 'EdgeColor', 'g', 'LineWidth', 3);
end
axis image xy off;
colorbar;
colormap('gray');
% title('Compute Subaperture Reference Positions');
hold off;
drawnow;

% exportgraphics(fig1, 'FLIR_Reference_Frame.pdf');

%% compute spot centers (centroid)

warning('off');

tic;

% pre-allocate for centroid position
xc_imgs_flir = NaN(length(xc_ref_flir), size(imgs_thresh, 3));
yc_imgs_flir = NaN(size(xc_imgs_flir));

% loop over frames
for ii = 1:size(xc_imgs_flir, 2)

    fprintf(1, ['Processing FLIR image ', num2str(ii), ' of ', num2str(size(imgs_thresh, 3)), '.\n']);

    % loop over subapertures
    for jj=1:size(xc_imgs_flir, 1)

        % define mask for the subaperture
        mask = zeros(size(imgs_thresh, [1,2]));
        row1_tmp = yc_ref_flir(jj)-subap_sz_flir/2;
        row2_tmp = yc_ref_flir(jj)+subap_sz_flir/2-1;
        col1_tmp = xc_ref_flir(jj)-subap_sz_flir/2;
        col2_tmp = xc_ref_flir(jj)+subap_sz_flir/2-1;

        if row1_tmp <= 0
            row1_tmp = 1;
        end
        if row1_tmp > 1536
            row1_tmp = 1536;
        end
        if row2_tmp <= 0
            row2_tmp = 1;
        end
        if row2_tmp > 1536
            row2_tmp = 1536;
        end

        if col1_tmp <= 0
            col1_tmp = 1;
        end
        if col1_tmp > 2048
            col1_tmp = 2048;
        end
        if col2_tmp <= 0
            col2_tmp = 1;
        end
        if col2_tmp > 2048
            col2_tmp = 2048;
        end
        
        

        % mask(yc_ref_flir(jj)-subap_sz_flir/2:yc_ref_flir(jj)+subap_sz_flir/2-1,...
        %     xc_ref_flir(jj)-subap_sz_flir/2:xc_ref_flir(jj)+subap_sz_flir/2-1) = 1;
        mask(row1_tmp:row2_tmp, col1_tmp:col2_tmp) = 1;

        % masked image
        tmp = imgs_thresh(:,:,ii) .* mask;

        % value and subscript of each non-zero pixel
        [rc, cc, vc] = find(tmp);

        if ~isempty(vc) % if subaperture has signal
            [~, idx_plot_sort] = sort(vc);

            % c_test = winter(length(idx_plot_sort));
            % 
            % cc_centroid = sum(cc.*vc) / sum(vc);
            % rc_centroid = sum(rc.*vc) / sum(vc);
            % 
            % fig1 = figure;
            % scatter3(cc(idx_plot_sort), rc(idx_plot_sort), vc(idx_plot_sort), 8, c_test, 'filled');
            % hold on;
            % h = plot3([cc_centroid, cc_centroid], [rc_centroid, rc_centroid], [0, ceil(max(vc))+2], 'k-', 'LineWidth', 2);
            % xlabel('Detector Column (pixel)', 'FontWeight', 'bold');
            % ylabel('Detector Row (pixel)', 'FontWeight', 'bold');
            % zlabel('Intensity', 'FontWeight', 'bold');
            % legend(h, {'Centroid'}, 'location', 'northeast');
            % ax = gca;
            % ax.FontWeight = 'bold';
            % hold off;
            % 
            % exportgraphics(fig1, 'Centroid_Computation.pdf');

            % compute centroid of all non-zero pixels
            xc_imgs_flir(jj,ii) = sum(cc.*vc) / sum(vc);
            yc_imgs_flir(jj,ii) = sum(rc.*vc) / sum(vc);

        else % if subaperture does not have signal

            xc_imgs_flir(jj,ii) = NaN;
            yc_imgs_flir(jj,ii) = NaN;

        end

    end

end
clear('tmp');

fprintf(1, 'Done. '); toc;

warning('on');

clear('rc', 'yc', 'vc', 'mask');

%%

figure;
plot(1:size(xc_imgs_flir, 2), xc_imgs_flir, 'LineWidth', 1);

figure;
plot(1:size(yc_imgs_flir, 2), yc_imgs_flir, 'LineWidth', 1);

%%

flag_mkVideo = true;
if flag_mkVideo
    v = VideoWriter('Frames_TrackedSpots', 'MPEG-4');
    v.Quality = 100;
    v.FrameRate = 10;
    open(v);
end

fig1 = figure;
hold on;
% for ii = size(imgs_thresh, 3)-400:size(imgs_thresh, 3)
for ii = 200:500
    imagesc(imgs_thresh(:,:,ii));
    axis image xy;
    colormap(gray);
    colorbar;
    hold on;
    plot(xc_imgs_flir(:,ii), yc_imgs_flir(:,ii), 'rx', 'LineWidth', 2, 'MarkerSize', 5);
    % xlim([0, 700]);
    % ylim([400, 1100]);
    title(['Frame ', num2str(ii, '%04d')]);

    drawnow;
    if flag_mkVideo
        writeVideo(v, getframe(fig1));
    end
    clf(fig1);
end
close(fig1);

if flag_mkVideo
    close(v);
end

%%

xc_imgs_flir_mean = mean(xc_imgs_flir(:,1:500), 2);
yc_imgs_flir_mean = mean(yc_imgs_flir(:,1:500), 2);

fig1 = figure;
histogram(xc_imgs_flir(:,1:500) - xc_imgs_flir_mean);
hold on;
histogram(yc_imgs_flir(:,1:500) - yc_imgs_flir_mean);
xlabel('Slope (Pixels)');
ylabel('Counts');
legend('x-axis', 'y-axis');
hold off;

% exportgraphics(fig1, 'Slope_Histogram.pdf');

%% save

filename_save = ['tracks_flir_dataset', num2str(dataset_flag, '%02d'), '.mat'];
save(filename_save, 'xc_imgs_flir', 'yc_imgs_flir', 'thresh_global', 'thresh_reference', 'flirFrame_ref');