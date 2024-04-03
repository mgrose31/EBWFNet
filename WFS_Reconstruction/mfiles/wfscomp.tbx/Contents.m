% Wfscomp Toolbox (WFSCOMP)
% Version 1.0.31 18-Jun-2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WFSCOMP TOOLS
%
% MZA Associates Corporation
% 6651 Centerville Business Pkwy Ste B
% Dayton OH 45459-2678
%
% Send bug reports to:  emagee@mza.com, 937-432-6560 ext. 243
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
%
% Type HELP 'function name' for more details on listed functions. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions:
% aowfsAccumFrame - accumulates single frame quantities into sums over 
%                   all frames of data set
% aowfsComputeStats - calculates results for the entire data set
% aowfsGetFrameInfo - applies aperture mask and calculates single frame 
%                     slope and irradiance values
% aowfsGetFrmstrs - prepares character strings for accessing data variable 
%                   names by frame number
% aowfsIrrVarSW - calculates sliding window averages of the irradiance 
%                 variance 
% aowfsGetProcInfo - sets up parameters and subaperture position data for 
%                    further data processing 
% aowfsSlpStrFn - locates correct center radius needed for conversion.
% ------------------------------------------------------------------------
% Functions not specific to aowfs, used in similiar form in other processing:
% cpaperstart - produces values of x_self slop struct function
% datetime  - puts a timestamp into a matlab date vector form
% getR0input - 
% makemask  - sets up an aperture mask vector from given parameters
% mknollv - 
% mkslpstfb1 - locates correct center radius
% mkwv - computes zernike coefficent to slopes matrix w.
% pzernikev - 
% pzrntonmv - computes mapping from an indexing of the Zernike polynomials.
% r0slpstrf - computes experimental and theoretical slope struct functions
% r0test - 
% r0tosc - analyze the r0 of the data contained in outz 
% ------------------------------------------------------------------------
% Data Structure Documentation (private):
% procInfo  - basic processing parameters and subaperture locations
% frameInfo  - slopes, irradiances, Zernikes, and slope structure 
%                 functions for each frame
% accum  - accumulated sums, means, and variances for frames of data set
% stats  - primary analysis results for data set
% ------------------------------------------------------------------------
% Files (private):
% aodueck1.mat  - test data file
% aodueck2.mat  - test data file
% aotoolscript.txt  - script for using aotool
% aotoscapctrs.mat  - subaperture center locations (from toscwfs.mat)
% aowfsINSTRUCTIONS.txt  - controls and executes aowfs data processing
% apermask.mat  - aperture mask vector made from data in toscwfs.mat
% toscwfs.mat  - fundamental parameters of aowfs apparatus 
% trialfun - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%