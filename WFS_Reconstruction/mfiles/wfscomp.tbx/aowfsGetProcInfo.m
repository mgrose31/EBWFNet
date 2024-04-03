function procInfo = aowfsGetProcInfo(nz);
% SYNTAX: 
% procInfo = aowfsGetProcInfo(nz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% aowfsGetProcInfo computes quantities necessary for processing the aowfs 
% data. It loads the locations of the subaperture centers 
% (aotoscsapctrs.mat) and the aperture mask (apermask.mat, a column vector 
% of the non-zero subaperture indices) from files which have been placed 
% in the working directory.  Other necessary parameters are specified 
% or calculated, and placed in the output data structure, procInfo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% nz [ ] =  The desired number of Zernike terms.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% procInfo [ ] =  A data structure containing basic processing parameters,
%                 documented in procInfo.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: aowfsGetProcInfo.m 3051 2010-10-01 20:33:26Z amoran $

%% BEGIN_CODE
 
load aotoscsapctrs
load apermask 
xsubctrs = ysubtosc;
ysubctrs = xsubtosc;
nsubs = 16;
xsubaps = 16;
ysubaps = 16;
nsubaps = nsubs^2;
spacing =  0.0475;  % subap width in meters
inradius = (1.4)*spacing; 
outradius = (8.6)* spacing;
xcenter = 0.0;
ycenter = 0.0;
xsubcor = xsubctrs - spacing/2;
ysubcor = ysubctrs - spacing/2;

% Note:  For r0 calculations, radii must include corners of subap or routines will
% drop the subap from calculations. (1.4) is chosen as slightly less than 2^(1/2), and
% (8.6) as slightly greater than 6*this, based on the aotool figure which shows inner radius
% as 1 subap diagonal and outer radius as 6 subap diagonals.  The exact numbers from the aotool
% figure are inner=0.083206 and outer=0.3875, vs 0.067175 and 0.40305 used by code.

% The following s structure is used below in calculating Zernike coefficients
clear s;
s.spacing = spacing;
s.xcenter = xcenter;
s.ycenter = ycenter;
s.outradius = outradius;
s.xsub = xsubcor; %THIS IS NOT THE SAME AS XSUB STORED IN PROCINFO BELOW--correct this 
s.ysub = ysubcor; %THIS IS NOT THE SAME AS YSUB STORED IN PROCINFO BELOW--correct this

countstorad = (1.0e-06/(16.0*2.3/3.0))/spacing;
ivtocn2ryt = (5.3056e-17)/(4 * 0.3090787800609279);  % parameter for cn2 calculation
framerate = 2300.0;

procInfo.ivtocn2ryt = ivtocn2ryt;
procInfo.countstorad = countstorad;
procInfo.framerate = framerate;
procInfo.lambda = 0.995e-06;
procInfo.nsubaps = nsubaps;
procInfo.apermask = apermask;

procInfo.xcenter = xcenter;
procInfo.ycenter = ycenter;
procInfo.spacing = spacing;
procInfo.inradius = inradius;
procInfo.outradius = outradius;
procInfo.xsub = xsubctrs;
procInfo.ysub = ysubctrs;

% calculate and store Zernike w matrix using s structure created above

if (nz ~= 0)
   w = mkwv(s,nz,6);
   winv = pinv(w);

   procInfo.w = w;
   procInfo.winv = winv;
end


