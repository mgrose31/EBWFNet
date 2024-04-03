clear; close all; clc;
addpath('..\..\Filtering');

% Code to match subapertures between event and frame cameras. For
% performing wavefront reconstruction with MZA's protocol, select valid
% subapertures from bottom to top an then from left to right.
% Author: Mitchell Grose, University of Dayton, 2023

%% set subaperture sizes

% subaperture sizes
subap_sz_ev = 54; % Prophesee VGA sensor
subap_sz_flir = 200; % FLIR camera

%% load event data

date_flag = '2023-07-30';
dataset_flag = 35;

fpath_ev = ['D:\Data\Dual_WFS\', date_flag, '_Eastwood\DualSync\Event\dataset', num2str(dataset_flag, '%02d')];
fpath_flir = ['D:\Data\Dual_WFS\', date_flag, '_Eastwood\DualSync\FLIR\dataset', num2str(dataset_flag, '%02d')];

% load event data
load(fullfile(fpath_ev, 'eventData_cd.mat'));

%% process event data

% add 1 to events due to MATLAB indexing
events.x = events.x + 1;
events.y = events.y + 1;

idx = events.ts >= events.trigger.ts(1) & events.ts <= events.trigger.ts(end);
events.x = events.x(idx);
events.y = events.y(idx);
events.p = events.p(idx);
events.ts = events.ts(idx);

% remove hot pixels
flag_removeHotPixels = true;
if flag_removeHotPixels
    tic; fprintf(1, 'Removing hot pixels... ');

    load('positions_hot_pixels.mat');

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

% inceptive event (IE) filter
flag_IE_filter = true;
if flag_IE_filter
    tic;

    twindow_IE_filter = 5000;

    fprintf(1, ['Event filtering with ', num2str(twindow_IE_filter), ' microsecond window... ']);

    [IE, IEm] = IE_filter(events, twindow_IE_filter);

%     IE = IE & (IEm > 1);

    events.x = events.x(IE);
    events.y = events.y(IE);
    events.p = events.p(IE);
    events.ts = events.ts(IE);
    events.IEm = IEm(IE);

    fprintf(1,'Done! '); toc;
else
    twindow_IE_filter = [];
end

height = events.height;
width = events.width;

% accumulate events to build an event reference frame
evFrame_ref = accumarray([events.y, events.x], 1, [height, width], @sum, 0);

%% load FLIR data and compute reference frame

dd_FLIR = dir(fullfile(fpath_flir, '*.tiff'));

numFramesLoad = length(dd_FLIR);

% load images and format
flirFrame_ref = zeros([1536, 2048]);

fprintf(1, 'reading FLIR images... ');
for ii=1:numFramesLoad
    % fliplr because 45 degree beamsplitter flips the image left-right
    flirFrame_ref = flirFrame_ref + fliplr(double(imread(fullfile(dd_FLIR(ii).folder, dd_FLIR(ii).name))));
end

flirFrame_ref = flirFrame_ref / numFramesLoad;

thresh_reference = 0.4;

flirFrame_ref = flirFrame_ref - max(flirFrame_ref(:)) .* thresh_reference;
flirFrame_ref(flirFrame_ref < 0) = 0;

%% select event spots

% plot FLIR reference image for selecting event spot locations
figure;
imagesc(flirFrame_ref);
axis image xy;
colormap(gray); colorbar;

ref_pos_ev_filename = ['ref_pos_ev_dataset', num2str(dataset_flag, '%02d'), '.mat'];
ref_pos_ev_filename_tmp = ['ref_pos_ev_dataset', num2str(dataset_flag, '%02d'), '_tmp.mat'];

if isfile(ref_pos_ev_filename_tmp)
    load(ref_pos_ev_filename_tmp)

    figure;
    imagesc(evFrame_ref);
    axis image xy;
    colormap(gray);
    colorbar;
    hold on;
    plot(xc_ref_ev_tmp, yc_ref_ev_tmp, 'rx', 'LineWidth', 2);
    hold off;

else
    fprintf(1, 'Select event spot centers!\n');

    figure;
    imagesc(evFrame_ref);
    axis image xy;
    colormap(gray);
    colorbar;

    [xc_ref_ev_tmp, yc_ref_ev_tmp] = ginput();
    save(ref_pos_ev_filename_tmp, 'xc_ref_ev_tmp', 'yc_ref_ev_tmp');
end

%% select FLIR spots

close all;

figure;
imagesc(evFrame_ref);
axis image xy;
colormap(gray);
colorbar;
hold on;
plot(xc_ref_ev_tmp, yc_ref_ev_tmp, 'rx', 'LineWidth', 2);
hold off;

ref_pos_flir_filename = ['ref_pos_flir_dataset', num2str(dataset_flag, '%02d'), '.mat'];
ref_pos_flir_filename_tmp = ['ref_pos_flir_dataset', num2str(dataset_flag, '%02d'), '_tmp.mat'];


if isfile(ref_pos_flir_filename_tmp)
    load(ref_pos_flir_filename_tmp);

    figure;
    imagesc(flirFrame_ref);
    axis image xy;
    colormap(gray);
    colorbar;
    hold on;
    plot(xc_ref_flir_tmp, yc_ref_flir_tmp, 'rx', 'LineWidth', 2);
    hold off;

else
    figure;
    imagesc(flirFrame_ref);
    axis image xy;
    colormap(gray);
    colorbar;

    [xc_ref_flir_tmp, yc_ref_flir_tmp] = ginput();

    save(ref_pos_flir_filename_tmp, 'xc_ref_flir_tmp', 'yc_ref_flir_tmp');
end

figure;
imagesc(flirFrame_ref);
axis image xy;
colormap(gray);
colorbar;
hold on;
plot(xc_ref_flir_tmp, yc_ref_flir_tmp, 'rx', 'LineWidth', 2);
hold off;

%% use template match to find center of event sub-apertures

% build meshgrid for template
[xx, yy] = meshgrid(-20:20, -20:20);
[theta, r] = cart2pol(xx, yy);
theta = rad2deg(theta);

% use a disk
% r(r>=9) = 0; r(r<5) = 0;
r(r>15) = 0;
h = double(logical(r));
h(21,21) = 1;


% size of the kernel
[M,N] = size(h);

% cross-correlate template and event frames
C = normxcorr2(h, evFrame_ref);
C = C(M/2:end-M/2, N/2:end-N/2); % crop image to original size

% figure;
% imagesc(C);
% axis image;
% colormap(turbo);
% colorbar;

% pre-allocate arrays to hold x,y position and template match correlation score
xc_imgs_ev = zeros(length(xc_ref_ev_tmp), 1);
yc_imgs_ev = zeros(size(xc_imgs_ev));
vc_imgs_ev = zeros(size(xc_imgs_ev));

% loop over the subapertures
for jj=1:length(yc_ref_ev_tmp)

    % set mask
    mask = zeros(size(evFrame_ref, [1,2]));
    mask(yc_ref_ev_tmp(jj)-subap_sz_ev/2:yc_ref_ev_tmp(jj)+subap_sz_ev/2-1, xc_ref_ev_tmp(jj)-subap_sz_ev/2:xc_ref_ev_tmp(jj)+subap_sz_ev/2-1) = 1;

    % apply mask to cross-correlation images
    C_tmp = C .* mask;

    [vc, I_tmp] = max(C_tmp(:));
    [rc, cc] = ind2sub(size(C_tmp), I_tmp);

    % yc_imgs_ev(jj,ii) = rc(1);
    % xc_imgs_ev(jj,ii) = cc(1);
    % vc_imgs_ev(jj,ii) = vc(1);

    % use a 3x3 neighborhood weighted average
    rc_arr = rc(1)-4:rc(1)+4;
    cc_arr = cc(1)-4:cc(1)+4;
    vc_arr = C_tmp(rc_arr, cc_arr);
    xc_imgs_ev(jj,1) = sum(repmat(cc_arr, [length(rc_arr),1]) .* (vc_arr ./ sum(vc_arr(:))), 'all');
    yc_imgs_ev(jj,1) = sum(repmat(rc_arr', [1,length(rc_arr)]) .* (vc_arr ./ sum(vc_arr(:))), 'all');
    vc_imgs_ev(jj,1) = mean(vc_arr(:));

end

xc_ref_ev = xc_imgs_ev;
yc_ref_ev = yc_imgs_ev;

figure;
imagesc(C);
axis image xy;
colormap(turbo);
colorbar;
hold on;
plot(xc_ref_ev, yc_ref_ev, 'kx', 'LineWidth', 2);
hold off;

figure;
imagesc(evFrame_ref);
axis image xy;
colormap(gray);
colorbar;
hold on;
plot(xc_ref_ev, yc_ref_ev, 'rx', 'LineWidth', 2);
hold off;

% save the center of the event subapertures
save(ref_pos_ev_filename, 'xc_ref_ev', 'yc_ref_ev');


%% compute spot centers (centroid)

% pre-allocate for centroid position
xc_ref_flir = NaN(length(xc_ref_flir_tmp), 1);
yc_ref_flir = NaN(size(xc_ref_flir_tmp));

% loop over subapertures
for jj=1:size(xc_ref_flir_tmp, 1)

    % define mask for the subaperture
    mask = zeros(size(flirFrame_ref, [1,2]));
    mask(yc_ref_flir_tmp(jj)-subap_sz_flir/2:yc_ref_flir_tmp(jj)+subap_sz_flir/2-1,...
        xc_ref_flir_tmp(jj)-subap_sz_flir/2:xc_ref_flir_tmp(jj)+subap_sz_flir/2-1) = 1;

    % masked image
    tmp = flirFrame_ref .* mask;

    % value and subscript of each non-zero pixel
    [rc, cc, vc] = find(tmp);

    if ~isempty(vc) % if subaperture has signal

        % compute centroid of all non-zero pixels
        xc_ref_flir(jj,1) = sum(cc.*vc) / sum(vc);
        yc_ref_flir(jj,1) = sum(rc.*vc) / sum(vc);

    else % if subaperture does not have signal

        xc_ref_flir(jj,1) = NaN;
        yc_ref_flir(jj,1) = NaN;

    end
end

figure;
imagesc(flirFrame_ref);
axis image xy;
colormap(gray);
colorbar;
hold on;
plot(xc_ref_flir, yc_ref_flir, 'rx', 'LineWidth', 2);
hold off;

% save the center of the frame (FLIR) subapertures
save(ref_pos_flir_filename, 'xc_ref_flir', 'yc_ref_flir');

%% transform from frame to events (TEST CODE)
% 
% x1 = xc_ref_flir;
% y1 = yc_ref_flir;
% x2 = xc_ref_ev;
% y2 = yc_ref_ev;
% 
% % transform FLIR frames to event frames
% tform_FLIR2ev = fitgeotform2d([x1, y1], [x2, y2], 'projective');
% 
% [x1T, y1T] = transformPointsForward(tform_FLIR2ev, x1, y1);
% 
% % transform FLIR tracks to event space
% % [xc_imgs_flirT, yc_imgs_flirT] = transformPointsForward(tform_FLIR2ev, xc_imgs_flir, yc_imgs_flir);
% img = imwarp(flirFrame_ref, tform_FLIR2ev, OutputView = imref2d(size(evFrame_ref)));
% 
% figure;
% imagesc(evFrame_ref);
% axis image xy;
% colormap(gray);
% colorbar;
% hold on;
% plot(xc_imgs_ev, yc_imgs_ev, 'rx', 'LineWidth', 2);
% plot(x1T, y1T, 'gx', 'LineWidth', 2);
% legend('Event', 'Flir Transformed');
% hold off;
% 
% 
% hf = figure;
% h1 = axes;
% imagesc(img);
% % xlim([50,250]);
% % ylim([50,250]);
% colormap(h1, 'gray');
% axis image xy;
% 
% h2 = axes;
% imagesc(evFrame_ref, 'AlphaData', evFrame_ref > 100);
% % xlim([50,250]);
% % ylim([50,250]);
% set(h2, 'color', 'none', 'visible', 'off');
% axis image xy;
% 
% % h2 = axes;
% % imagesc(B_pos, 'AlphaData', B_pos > 0.5);
% % set(h2, 'color', 'none', 'visible', 'off');
% % colormap(h2, [1, 0, 0]);
% % 
% % h3 = axes;
% % imagesc(B_neg, 'AlphaData', B_neg < -0.5);
% % set(h3, 'color', 'none', 'visible', 'off');
% % colormap(h3, [0, 0, 1]);
% % linkaxes([h1, h2, h3]);