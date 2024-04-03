function frameInfo = aowfsGetFrameInfo2(irigin,procInfo,nz,flagR0slp,status);
% SYNTAX: 
% frameInfo = aowfsGetFrameInfo2(irigin,procInfo,nz,flagR0slp,status);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:  
% aowfsGetFrameInfo applies the aperture mask to the spot position 
% (ie. slope) and  intensity (ie. irradiance) data for each frame.  Because 
% these quantities are stored in  separately named variables for each 
% frame, some string manipulation is necessary to access them,  which 
% requires the function aowfsGetFrmstrs.  Slope and irradiance quantities 
% are calculated for each frame and placed in the data structure frameInfo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% nz [ ] =  The desired number of Zernike terms.
% flagR0slp [ ] =  A flag for the calculation of r0 from slope structure 
%                  functions. flagR0slp=1 enables this calculation, 
%                  flagR0slp=0 omits it. The slope structure functions are 
%                  calculated for each frame by aowfsGetFrameInfo, but the 
%                  r0 calculation is only done for the entire data set by 
%                  aowfsComputeStats.
% filename [ ] =  A string specifying complete file name for the data set.
% procInfo [ ] =  The data structure containing basic processing parameters 
%                 output by aowfsGetProcInfo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% frameInfo [ ] =  An array of data structures containing slope, 
%                  irradiance, and cn2 quantities calculated for each 
%                  frame, documented in frameInfo.m
%
% NOTE: frmstrs is a cell array of strings utilized to process quantites 
%       on a per frame basis, using eval. frmstrs is the output of the 
%       function aowfsGetFrmstrs, which MUST BE MODIFIED IF THE NAMING 
%       CONVENTION FOR THE VARIABLES CONTAINING THE DATA IS CHANGED.  The 
%       eval statements in   aowfsGetFrameInfo below may also need 
%       modification in such circumstances.
%
% There are separately named variables for each frame: "xspotpos0000003", 
% etc. frmstrs{ii} contains the seven final numeric characters, which 
% specify the frame number, ii. Because frame 1 data is specified by 
% '0000000',  frame 224 is specified by '0000223', etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: aowfsGetFrameInfo2.m 3051 2010-10-01 20:33:26Z amoran $

%% BEGIN_CODE
 
if (exist('status','var'))
   if (status ~= 0)
      tic;
      tic;
   end;
end;
load(['int',irigin]);
load(['slp',irigin]);
numframes = nframes;
irig0 = convirig(irig); % initial time in days.
if (~exist('framerate','var'))
   framerate = procInfo.framerate;
end;
irigdt = (60*60*24)*(1.0/framerate); % dt in days.
if (exist('MicronsPerSlopeBit','var'))
   procInfo.countstorad = (MicronsPerSlopeBit*1.0e-06)/procInfo.spacing;
end;
frmstrs = aowfsGetFrmstrs(numframes);
nsub = procInfo.nsubaps;
countstorad = procInfo.countstorad;
ivtocn2ryt = procInfo.ivtocn2ryt;  % parameter for cn2 calculation
apermask = procInfo.apermask;
if (exist('status','var'))
   if (status ~= 0)
      'prepare', toc,
   end;
end;

if (exist('status','var'))
   if (status ~= 0)
      tic;
   end;
end;
for ii=1:1:numframes
   if (exist('status','var'))
      if (status > 1)
         tic;
      end;
   end;
   eval(['xpos=xsp',frmstrs{ii},';']);
   eval(['ypos=ysp',frmstrs{ii},';']);
   eval(['inten=int',frmstrs{ii},';']);
   clear(['xsp',frmstrs{ii}], ['ysp',frmstrs{ii}],['int',frmstrs{ii}]);

   frameInfo(ii).time = irig0 + (ii-1)*irigdt;

   slopex = reshape(xpos,nsub,1)*countstorad;
   slopey = reshape(ypos,nsub,1)*countstorad;
   irr = reshape(inten,nsub,1);
   slopex = slopex(apermask);
   slopey = slopey(apermask);
   irr = irr(apermask);
   frameInfo(ii).slopex = slopex;
   frameInfo(ii).slopey = slopey;
   frameInfo(ii).slopexBar = mean(frameInfo(ii).slopex);
   frameInfo(ii).slopeyBar = mean(frameInfo(ii).slopey);

%
%  compute irradiance quantities:
%
   irrBar = mean(irr);
   irr2Bar = mean(irr.^2);
   irrVar = (std(irr))^2;
   irrFNorm = irr/irrBar;
   irrFNVar = (std(irrFNorm))^2;
   cn2FN = ivtocn2ryt * irrFNVar;

   frameInfo(ii).irr = irr;
   frameInfo(ii).irrBar = irrBar;
   frameInfo(ii).irr2Bar = irr2Bar;
   frameInfo(ii).irrVar = irrVar;
   frameInfo(ii).irrFNVar = irrFNVar;
   frameInfo(ii).cn2FN = cn2FN;

%
%  compute Zernikes
%
   if (nz ~= 0)
       slopes = [slopex' slopey'];
       frameInfo(ii).zern = procInfo.winv*slopes.';
   end

%
%  compute tilt structure functions
%
   if (flagR0slp ~= 0)
       [dxxx,dxxy,dyyx,dyyy,dxxxy,dyyxy]=mkslpstfb1(slopex,slopey,...
                                                    procInfo.xsub,procInfo.ysub,...
                                                    procInfo.xcenter,procInfo.ycenter,...
                                                    procInfo.spacing,procInfo.inradius);
       frameInfo(ii).dxxx = dxxx.';
       frameInfo(ii).dxxy = dxxy.';
       frameInfo(ii).dyyx = dyyx.';
       frameInfo(ii).dyyy = dyyy.';
       frameInfo(ii).dxxxy = dxxxy.';
       frameInfo(ii).dyyxy = dyyxy.';
   end

   if ((ii == 1) & (numframes > 1))
      frameInfo(numframes) = frameInfo(1);
   end;

   if (exist('status','var'))
      if (status > 1)
         ii, toc,
      end;
   end;
end
if (exist('status','var'))
   if (status ~= 0)
      'loop', toc,
   end;
end;
if (exist('status','var'))
   if (status ~= 0)
      'completion', toc,
   end;
end;
