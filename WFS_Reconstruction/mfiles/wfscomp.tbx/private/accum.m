% accum:  A data structure containing accumulated slope and irradiance quantities summed over
%         all frames of the data set.  The output of function aowfsAccumFrame.
%    Note:  "Sum" indicates a sum over all frames processed, whereas
%           "Bar" indicates a mean value over a single frame of data
%
% accum.nframe:  The number of frames processed in the data set
%
% accum.slopexSum:  Sum of x slopes in each subaperture over all frames
% accum.slopeySum:  Sum of y slopes in each subaperture over all frames
% accum.slopexBarSum:  Sum of the single-frame mean values of x slopes
% accum.slopeyBarSum:  Sum of the single-frame mean values of y slopes
%
% accum.slopexSqSum:  Sum of squares of x slopes in each subaperture
% accum.slopeySqSum:  Sum of squares of y slopes in each subaperture
% accum.slopexBarSqSum:  Sum of squares of single-frame mean values of x slopes
% accum.slopeyBarSqSum:  Sum of squares of single-frame mean values of y slopes
%
% accum.irrSum:  Sum of irradiance in each subaperture over all frames
% accum.irrSqSum:  Sum of square of irradiance in each subaperture over all frames
% accum.irrBarSum:  Sum of means of single-frame irradiances
% accum.irrSqBarSum:  Sum of means of single-frame irradiance-squared values
% accum.irrVarSum:  Sum of variances of single-frame irradiances
%
% accum.cn2FNSum:  Sum of cn2 values calculated for each frame from 
%                  frame-normalized irradiances.
%
% accum.zernSum:  Sum of Zernike coefficients over all frames
% accum.zernSqSum:  Sum of squares of Zernike coefficients over all frames
%
% accum.dxxxSum:  Sum over all frames of x slope structure function in x direction 
% accum.dxxySum:  Sum over all frames of x slope structure function in y direction 
% accum.dyyxSum:  Sum over all frames of y slope structure function in x direction 
% accum.dyyySum:  Sum over all frames of y slope structure function in y direction 
% accum.dxxxySum:  Sum over all frames of x slope structure function in xy (diagonal) direction 
% accum.dyyxySum:  Sum over all frames of y slope structure function in xy (diagonal) direction 
%