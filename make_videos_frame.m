clear; close all; clc;

% pixel pitch of each sensor
pitch_EV = 15e-6; % Prophesee VGA, 15 micrometer pitch
pitch_FLIR = 3.45e-6; % FLIR, 3.45 micrometer pitch

date_flag = '2023-07-30';

thresh_global = 0.25;

%% load FLIR dataset

% for dataset_flag = 1:35
for dataset_flag = 4

    tic; fprintf(1, ['Processing ', num2str(dataset_flag), '... ']);

    fpath_flir = ['D:\Data\Dual_WFS\', date_flag, '_Eastwood\DualSync\FLIR\dataset', num2str(dataset_flag, '%02d')];

    dd_FLIR = dir(fullfile(fpath_flir, '*.tiff'));

    numFramesLoad = length(dd_FLIR);

    % load images and format
    imgs = zeros([1536, 2048, numFramesLoad], 'uint8');

    fprintf(1, 'reading FLIR images... ');
    for ii=1:numFramesLoad
        % fliplr because 45 degree beamsplitter flips the image left-right
        % imgs(:,:,ii) = fliplr(double(imread(fullfile(dd_FLIR(ii).folder, dd_FLIR(ii).name))));
        imgs(:,:,ii) = fliplr(imread(fullfile(dd_FLIR(ii).folder, dd_FLIR(ii).name)));
    end

    fprintf(1, 'writing raw video... ');

    v = VideoWriter(['dataset', num2str(dataset_flag, '%02d'), '_FLIR_Raw'], 'MPEG-4');
    v.Quality = 75;
    v.FrameRate = 10;

    open(v);

    fig1 = figure;
    hold on;
    axis image;
    colormap('gray');
    % for ii = 1:size(imgs, 3)
    for ii = 200:500
        img_tmp = imgs(:, :, ii);
        imagesc(img_tmp);
        colorbar; clim([min(img_tmp(:)), max(img_tmp(:))]);
        title(['Frame ', num2str(ii, '%04d')]);
        writeVideo(v, getframe(fig1));
        clf(fig1);
    end
    close(fig1);
    close(v);

    

    % apply global threshold to images
    fprintf(1, 'filtering images with global threshold... ');

    for ii=1:size(imgs, 3)
        img_tmp = imgs(:,:,ii);
        img_tmp = img_tmp - max(img_tmp, [], 'all').*thresh_global;
        img_tmp(img_tmp<0) = 0;
        imgs(:,:,ii) = img_tmp;        
    end

    fprintf(1, 'writing threshold video... ');

    v = VideoWriter(['dataset', num2str(dataset_flag, '%02d'), '_FLIR_Thresh'], 'MPEG-4');
    v.Quality = 75;
    v.FrameRate = 10;

    open(v);

    fig1 = figure;
    hold on;
    axis image;
    colormap('gray');
    for ii = 200:500
    % for ii = 1:size(imgs, 3)
    % for ii = 1500:size(imgs, 3)
        img_tmp = imgs(:, :, ii);
        imagesc(img_tmp);
        colorbar; clim([min(img_tmp(:)), max(img_tmp(:))]);
        title(['Frame ', num2str(ii, '%04d')]);
        writeVideo(v, getframe(fig1));
        clf(fig1);
    end
    close(fig1);
    close(v);


    clear('imgs', 'img_tmp');

    fprintf(1, 'Done. \n'); toc;

end

%%
% %% load FLIR dataset
% 
% % subaperture sizes
% subap_sz_flir = 200;
% 
% fpath_flir = ['D:\Data\Dual_WFS\', date_flag, '_Eastwood\DualSync\FLIR\dataset', num2str(dataset_flag, '%02d')];
% 
% dd_FLIR = dir(fullfile(fpath_flir, '*.tiff'));
% 
% if numFramesLoad == -1
%     numFramesLoad = length(dd_FLIR);
% end
% 
% % load images and format
% % imgs = NaN(1536, 2048, length(numFramesLoad));
% imgs = zeros([1536, 2048, numFramesLoad], 'uint8');
% 
% tic; fprintf(1, 'Reading FLIR images... \n');
% for ii=1:numFramesLoad
%     if rem(ii, 10) == 0
%         fprintf(1, ['Frame ', num2str(ii), ' of ', num2str(numFramesLoad), '.\n']);
%     end
%     % fliplr because 45 degree beamsplitter flips the image left-right
%     % imgs(:,:,ii) = fliplr(double(imread(fullfile(dd_FLIR(ii).folder, dd_FLIR(ii).name))));
%     imgs(:,:,ii) = fliplr(imread(fullfile(dd_FLIR(ii).folder, dd_FLIR(ii).name)));
% end
% fprintf(1, 'Done. '); toc;
% 
% % imgs = double(imgs);
% 
% %% compute reference image
% 
% thresh_reference = 0.5;
% 
% % compute reference image
% img_tavg = mean(imgs, 3);
% img_ref_flir = img_tavg - max(img_tavg(:)) .* thresh_reference;
% img_ref_flir(img_ref_flir < 0) = 0;
% 
% clear('img_tavg');
% 
% % plot reference image
% figure;
% imagesc(img_ref_flir);
% colormap('gray');
% colorbar;
% axis image; axis off;
% title('Intensity Sensor Reference Image');
% 
% %% make a video
% 
% v = VideoWriter('TestFLIR3', 'MPEG-4');
% v.Quality = 75;
% v.FrameRate = 20;
% 
% open(v);
% 
% fig1 = figure;
% hold on;
% axis image;
% colormap('gray');
% % for ii = 1:size(imgs, 3)
% for ii = 500:700
%     img_tmp = imgs(:, :, ii);
%     imagesc(img_tmp);
%     colorbar; clim([min(img_tmp(:)), max(img_tmp(:))]);
%     title(['Frame ', num2str(ii, '%04d')]);
%     writeVideo(v, getframe(fig1));
%     clf(fig1);
% end
% close(fig1);
% close(v);
