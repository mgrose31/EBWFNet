% frameInfo:  an array of data structures frameInfo(ii), the ouput of function aowfsGetFrameInfo,
%     containing the following fields for each processed frame (ii):
%
% frameInfo(ii).slopex:  subaperture spot slopes in the x direction--a column vector
% frameInfo(ii).slopey:  subaperture spot slopes in the y direction--a column vector
% frameInfo(ii).slopexBar:  mean of slopex for this frame
% frameInfo(ii).slopeyBar:  mean of slopey for this frame
% frameInfo(ii).irr:  subaperture irradiances--a column vector 
% frameInfo(ii).irrBar:  mean of subaperture irradiances for this frame 
% frameInfo(ii).irrVar:  variance of subaperture irradiances for this frame 
% frameInfo(ii).irrFNVar:  frame normalized variance of subaperture irradiances 
% frameInfo(ii).cn2FN:  cn2 calculated from frame normalized variance of subaperture irradiances 
% frameInfo(ii).zern:  matrix of Zernike coefficients for number of terms specified in nz parameter
%
% The following fields contain slope structure functions calculated over each frame from slope values.
%    In aowfsComputeStats these are used to estimate r0, if flagR0slp=1.
%    If flagR0slp=0, these fields are not computed and are set = 0.
%    Each of these is an average difference squared of x or y slope values, with the differences taken
%    between slope values for subapertures separated in the x, y, or xy directions.  Each structure function 
%    itself varies with subaperture separation, and is expressed as a vector of values for various subaperture
%    separations of 1, 2, ..., n subaperture widths.  Six different structure functions are obtained:
%
% frameInfo(ii).dxxx:  x slope structure function in x direction 
% frameInfo(ii).dxxy:  x slope structure function in y direction 
% frameInfo(ii).dyyx:  y slope structure function in x direction 
% frameInfo(ii).dyyy:  y slope structure function in y direction 
% frameInfo(ii).dxxxy:  x slope structure function in xy (diagonal) direction
% frameInfo(ii).dyyxy:  y slope structure function in xy (diagonal) direction 
%
