% stats:  A data structure containing principal analysis results for an entire data set, 
%         the output of function aowfsComputeStats.  It has the following fields:
%    Note: "Bar" indicates a mean value over a single frame of data, whereas
%          "Avg" indicates an average over all frames of the data set.
%
% stats.numframes:  The number of frames processed in the data set.
%
% stats.slopexAvg:  x slopes in each subaperture, averaged over all frames 
% stats.slopexVar:  variance over all frames of x slopes in each subaperture
% stats.slopeyAvg:  y slopes in each subaperture, averaged over all frames
% stats.slopeyVar:  variance over all frames of y slopes in each subaperture
%
% stats.slopexBarAvg:  mean of x slopes for each frame, averaged over all frames
% stats.slopexBarVar:  variance over all frames of mean of x slopes for each frame
% stats.slopeyBarAvg:  mean of y slopes for each frame, averaged over all frames
% stats.slopeyBarVar:  variance over all frames of mean of y slopes for each frame
%
% stats.irrAvgVar:  average irradiance variance for entire data set (not frame-normalized)
% stats.cn2Avg:  cn2 calculated from irrAvgVar for entire data set (not frame-normalized)
% stats.cn2FNAvg:  average over all frames of cn2 calculated per frame from frame-normalized irradiance
%
% stats.zernAvg:  Zernike coefficients averaged over all frames
% stats.zernVar:  variance over all frames of the Zernike coefficients
% stats.r0Zern:  r0 calculated from Zernike spectra
%
% stats.sfnXslpXdir:  x slope structure function in x direction
% stats.sfnXslpYdir:  x slope structure function in y direction
% stats.sfnYslpXdir:  y slope structure function in x direction
% stats.sfnYslpYdir:  y slope structure function in y direction
% stats.sfnXslpXYdir:  x slope structure function in xy direction
% stats.sfnYslpXYdir:  y slope structure function in xy direction
%
% stats.r0XslpXdir: r0 computed from x slope structure function in x direction
% stats.r0XslpYdir: r0 computed from x slope structure function in y direction
% stats.r0YslpXdir: r0 computed from y slope structure function in x direction
% stats.r0YslpYdir: r0 computed from y slope structure function in y direction
% stats.r0XslpXYdir:  r0 computed from x slope structure function in xy direction
% stats.r0YslpXYdir:  r0 computed from y slope structure function in xy direction
