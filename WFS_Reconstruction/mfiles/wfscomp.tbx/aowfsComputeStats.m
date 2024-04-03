function stats = aowfsComputeStats(accum,frameInfo,procInfo,nz,flagR0z,flagR0slp);
% SYNTAX: 
% stats = aowfsComputeStats(accum,frameInfo,procInfo,nz,flagR0z,flagR0slp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:  
% aowfsComputeStats uses the accumulated slope and irradiance data from a 
% wfs data set, along with Zernike coefficients and slope structure 
% functions, to compute various quantities of interest, including r0 and 
% cn2 estimates.  It produces the principal analysis results for an entire
% data set. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% accum [ ] =  The data structure containing accumulated slope and 
%              irradiance quantities output by aowfsAccumFrame.
% frameInfo [ ] =  The data structure containing slope and irradiance 
%                  quantities per frame, output by aowfsGetFrameInfo.
% procInfo [ ] =  The data structure containing basic processing parameters 
%                 output by aowfsGetProcInfo.
% nz [ ] =  The desired number of Zernike terms.
% flagR0z [ ] =  A flag for the calculation of R0 from Zernike coefficients.
% flagR0slp [ ] =  A flag for the calculation of r0 from slope structure 
%                  functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% stats [ ] =  A data structure containing principal analysis results for an
%              entire data set, documented in stats.m
% Note:  "Bar" refers to a mean value taken over one frame, "Avg" to an
% average over all frames of the data set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: aowfsComputeStats.m 3051 2010-10-01 20:33:26Z amoran $

%% BEGIN_CODE

nfr = accum.nframe;
stats.numframes = nfr;

stats.slopexAvg = accum.slopexSum/nfr;
stats.slopexVar = accum.slopexSqSum/nfr - stats.slopexAvg.^2;
stats.slopeyAvg = accum.slopeySum/nfr;
stats.slopeyVar = accum.slopeySqSum/nfr - stats.slopeyAvg.^2;

stats.slopexBarAvg = accum.slopexBarSum/nfr;
stats.slopexBarVar = accum.slopexBarSqSum/nfr - stats.slopexBarAvg.^2;
stats.slopeyBarAvg = accum.slopeyBarSum/nfr;
stats.slopeyBarVar = accum.slopeyBarSqSum/nfr - stats.slopeyBarAvg.^2;

% irradiance Avg Variance (not frame-normalized):

irrBarAvgSq = (accum.irrBarSum/nfr)^2;
irrSqBarAvg = accum.irrSqBarSum/nfr;
irrVarAvg = accum.irrVarSum/nfr;
irrAvgVar = (irrVarAvg + irrSqBarAvg - irrBarAvgSq)/irrBarAvgSq;
stats.irrAvgVar = irrAvgVar; 
stats.cn2Avg = procInfo.ivtocn2ryt * irrAvgVar;
stats.cn2FNAvg = accum.cn2FNSum/nfr;
if (nz~=0)
   stats.zernAvg = accum.zernSum/nfr;
   stats.zernVar = accum.zernSqSum/nfr - stats.zernAvg.^2;
   if (flagR0z~=0)
      [outz]=getR0input(frameInfo,stats,procInfo.lambda);
      [rout,sig,ssout]=r0tosc(outz,procInfo);
      stats.r0Zern = rout;
   end
end

if (flagR0slp~=0)
   if (~exist('outz','var'))
      [outz]=getR0input(frameInfo);
   end
   [r0XX,r0XY,r0XXY,r0YX,r0YY,r0YXY,dxxx,dxxy,dyyx,dyyy,dxxxy,dyyxy,dxxxth,dxxyth,dxxxyth]=r0slpstrf(outz,procInfo);

   stats.sfnXslpXdir = accum.dxxxSum/nfr;
   stats.sfnXslpYdir = accum.dxxySum/nfr;
   stats.sfnYslpXdir = accum.dyyxSum/nfr;
   stats.sfnYslpYdir = accum.dyyySum/nfr;
   stats.sfnXslpXYdir = accum.dxxxySum/nfr;
   stats.sfnYslpXYdir = accum.dyyxySum/nfr;

   stats.r0XslpXdir = r0XX;
   stats.r0XslpYdir = r0XY;
   stats.r0YslpXdir = r0YX;
   stats.r0YslpYdir = r0YY;
   stats.r0XslpXYdir = r0XXY;
   stats.r0YslpXYdir = r0YXY;

end

