clear; close all; clc;

% netH5 = load('./trained_network_depth4_slopeMSE_cartesian_fixedMaxPool.mat', 'net');
netH5 = load('trained_network_depth4_slopeMSE_cartesian_smallNetworkTest.mat', 'net');
load('.\Formatted_Data\formatted_dataset04.mat', 'subImgs', 'xx', 'yy');

net = netH5.net;
Layers = net.Layers;

layer = Layers(2).Name;
channels = 1:10;

X = subImgs(:,:,:,1,1);
X(:,:,end+1) = xx;
X(:,:,end+1) = yy;

X = dlarray(X, "SSC");

%%

% Y1 = forward(net, X, 'Outputs', layer);
Y1 = forward(net, X, 'Outputs', Layers(4).Name);
Y1 = extractdata(Y1);

figure('Position', [100, 100, 1000, 400]);
tiledlayout(3,5);

for ii = 1:15
    nexttile;
    imagesc(Y1(:,:,ii));
    axis image;
    colorbar;
    title(num2str(ii));
end