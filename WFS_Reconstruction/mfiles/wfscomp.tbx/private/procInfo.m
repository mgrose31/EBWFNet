% procInfo: a data structure containing basic processing parameters for wfs processing,
%     the output of function aowfsGetProcInfo, with the following fields:
%
% procInfo.ivtocn2ryt:   A parameter used in calculating cn2 from irradiance variances.
% procInfo.countstorad:  Conversion factor from slopes in meters at image plane to radians of tilt in object space.
% procInfo.nsubaps:  Total number of subapertures.  
% procInfo.apermask:  Aperture mask vector--a column vector of the non-zero subaperture indices.
% procInfo.xcenter:  x coordinate of the aperture center in meters.
% procInfo.ycenter:  y coordinate of the aperture center in meters.
% procInfo.spacing:  Spacing between subaperture centers in meters.
% procInfo.inradius:  Inner radius of the annular aperture in meters.
% procInfo.outradius:  Outer radius of the annular aperture in meters.
% procInfo.xsub:  x coordinates of the subaperture centers in meters--a column vector.
% procInfo.ysub:  y coordinates of the subaperture centers in meters--a column vector.
% procInfo.w:  The w matrix, for calculation of Zernike coefficients.
% procInfo.winv:  Inverse of the w matrix, for calculation of Zernike coefficients.
