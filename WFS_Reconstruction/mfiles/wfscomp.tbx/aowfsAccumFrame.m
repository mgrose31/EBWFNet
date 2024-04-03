function accum = aowfsAccumFrame(FrameInfo, first, last)
% SYNTAX:
% accum = aowfsAccumFrame(FrameInfo, first, last)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% aowfsAccumFrame sums data over all the frames of a data set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% FrameInfo [ ] =  The output data structure from aowfsGetFrameInfo, which
%                  is an array of structures indexed by frame number.  For 
%                  each frame, it contains slope and irradiance quantities 
%                  calculated for that frame.  The individual frames are 
%                  accessed by the assignment:  frameInfo=FrameInfo(ii); 
% first [ ] = First frame to be accumulated.
% last [ ] = Last frame to be accumulated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% accum [ ] =  A data structure containing accumulated slope and irradiance 
%              quantities summed over all frames of the data set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: aowfsAccumFrame.m 3051 2010-10-01 20:33:26Z amoran $

%% BEGIN_CODE

% Set up and zero fields:

if (missoa(FrameInfo))
   numframes = size(FrameInfo.time,2);
elseif (misaos(FrameInfo(1)))
   numframes = length(FrameInfo);
else
   accum = [];
   return;
end;
%
if (~exist('first','var'))
   first = 1;
end;
if (isempty(first))
   first = 1;
end;
%
if (~exist('last','var'))
   last = numframes;
end;
if (isempty(last))
   last = numframes;;
end;
%
if (isfield(FrameInfo(1), 'zern'))
   nz = 1;
end;
if (isfield(FrameInfo(1), 'dxxx'))
   flagR0slp = 1;
end;
%
if (~missoa(FrameInfo))
   accum.nframe = 0;
   accum.slopexSum = 0*FrameInfo(1).slopex;
   accum.slopeySum = accum.slopexSum;
   accum.slopexBarSum = 0;
   accum.slopeyBarSum = 0;
   accum.slopexSqSum = accum.slopexSum;
   accum.slopeySqSum = accum.slopexSum;
   accum.slopexBarSqSum = 0;
   accum.slopeyBarSqSum = 0;
   accum.irrSum = accum.slopexSum;
   accum.irrSqSum = accum.slopexSum; 
   accum.irrBarSum = 0;
   accum.irrSqBarSum = 0;
   accum.irrVarSum = 0;
   accum.cn2FNSum = 0;

   if (nz ~= 0)
      accum.zernSum = 0*FrameInfo(1).zern;
      accum.zernSqSum = accum.zernSum;
   end
   if (flagR0slp ~= 0)
      accum.dxxxSum = 0*FrameInfo(1).dxxx;
      accum.dxxySum = 0*FrameInfo(1).dxxy;
      accum.dyyxSum = 0*FrameInfo(1).dyyx;
      accum.dyyySum = 0*FrameInfo(1).dyyy;
      accum.dxxxySum = 0*FrameInfo(1).dxxxy;
      accum.dyyxySum = 0*FrameInfo(1).dyyxy;
   end
   
   %
   % Accumulate over all frames of data in the data set:
   %
   numframes = length(FrameInfo);
   for ii=first:1:last
      frameInfo=FrameInfo(ii);
      accum.nframe = accum.nframe + 1;
      accum.slopexSum = accum.slopexSum + frameInfo.slopex;
      accum.slopeySum = accum.slopeySum + frameInfo.slopey;
      accum.slopexBarSum = accum.slopexBarSum + frameInfo.slopexBar;
      accum.slopeyBarSum = accum.slopeyBarSum + frameInfo.slopeyBar;
      accum.slopexSqSum = accum.slopexSqSum + frameInfo.slopex.^2;
      accum.slopeySqSum = accum.slopeySqSum + frameInfo.slopey.^2;
      accum.slopexBarSqSum = accum.slopexBarSqSum + frameInfo.slopexBar.^2;
      accum.slopeyBarSqSum = accum.slopeyBarSqSum + frameInfo.slopeyBar.^2;
      accum.irrSum = accum.irrSum + frameInfo.irr;
      accum.irrSqSum = accum.irrSqSum + frameInfo.irr.^2;
      accum.irrBarSum = accum.irrBarSum + frameInfo.irrBar;
      accum.irrSqBarSum = accum.irrSqBarSum + mean(frameInfo.irr.^2);
      accum.irrVarSum = accum.irrVarSum + frameInfo.irrVar;
      accum.cn2FNSum = accum.cn2FNSum + frameInfo.cn2FN;
      if (nz ~= 0)
         accum.zernSum = accum.zernSum + frameInfo.zern;
         accum.zernSqSum = accum.zernSqSum + frameInfo.zern.^2;
      end
      if (flagR0slp ~=0)
         accum.dxxxSum = accum.dxxxSum + frameInfo.dxxx;
         accum.dxxySum = accum.dxxySum + frameInfo.dxxy;
         accum.dyyxSum = accum.dyyxSum + frameInfo.dyyx;
         accum.dyyySum = accum.dyyySum + frameInfo.dyyy;
         accum.dxxxySum = accum.dxxxySum + frameInfo.dxxxy;
         accum.dyyxySum = accum.dyyxySum + frameInfo.dyyxy;
      end
   end;
else
   accum.nframe = last - first + 1;
   accum.slopexSum = sum(FrameInfo.slopex(:,first:last).').';
   accum.slopeySum = sum(FrameInfo.slopey(:,first:last).').';
   accum.slopexBarSum = sum(FrameInfo.slopexBar(:,first:last).').';
   accum.slopeyBarSum = sum(FrameInfo.slopeyBar(:,first:last).').';
   accum.slopexSqSum = sum((FrameInfo.slopex(:,first:last).').^2).';
   accum.slopeySqSum = sum((FrameInfo.slopey(:,first:last).').^2).';
   accum.slopexBarSqSum = sum((FrameInfo.slopexBar(:,first:last).').^2).';
   accum.slopeyBarSqSum = sum((FrameInfo.slopeyBar(:,first:last).').^2).';
   accum.irrSum = sum(FrameInfo.irr(:,first:last).').';
   accum.irrSqSum = sum((FrameInfo.irr(:,first:last).').^2).';
   accum.irrBarSum = sum(FrameInfo.irrBar(:,first:last).').';
   accum.irrSqBarSum = sum(mean(FrameInfo.irr(:,first:last).^2).').';
   accum.irrVarSum = sum(FrameInfo.irrVar(:,first:last).').';
   accum.cn2FNSum = sum(FrameInfo.cn2FN(:,first:last).').';
   if (nz ~= 0)
      accum.zernSum = sum(FrameInfo.zern(:,first:last).').';
      accum.zernSqSum = sum((FrameInfo.zern(:,first:last).^2).').';
   end
   if (flagR0slp ~=0)
      accum.dxxxSum = sum(FrameInfo.dxxx(:,first:last).').';
      accum.dxxySum = sum(FrameInfo.dxxy(:,first:last).').';
      accum.dyyxSum = sum(FrameInfo.dyyx(:,first:last).').';
      accum.dyyySum = sum(FrameInfo.dyyy(:,first:last).').';
      accum.dxxxySum = sum(FrameInfo.dxxxy(:,first:last).').';
      accum.dyyxySum = sum(FrameInfo.dyyxy(:,first:last).').';
   end;
end
