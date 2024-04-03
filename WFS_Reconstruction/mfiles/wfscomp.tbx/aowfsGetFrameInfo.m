function frameInfo = aowfsGetFrameInfo(nz,flagR0slp,irigin,procInfo);
% SYNTAX:
% frameInfo = aowfsGetFrameInfo(nz,flagR0slp,irigin,procInfo);
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
% nz [ ] = The desired number of Zernike terms.
% flagR0slp [ ] =  A flag for the calculation of r0 from slope structure 
%                  functions. flagR0slp=1 enables this calculation, 
%                  flagR0slp=0 omits it. The slope structure functions are 
%                  calculated for each frame by aowfsGetFrameInfo, but the 
%                  r0 calculation is only done for the entire data set by 
%                  aowfsComputeStats.
% filename [ ] =  A string specifying the complete file name for the data set.
% procInfo [ ] =  The data structure containing basic processing parameters 
%                 output by aowfsGetProcInfo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% frameInfo [ ] =  An array of data structures containing slope, 
%                  irradiance, and cn2 quantities calculated for each 
%                  frame, documented in frameInfo.m
% NOTE: frmstrs is a cell array of strings utilized to process quantites on 
%       a per frame basis, using eval. frmstrs is the output of the 
%       function aowfsGetFrmstrs, which MUST BE MODIFIED IF THE NAMING
%       CONVENTION FOR THE VARIABLES CONTAINING THE DATA IS CHANGED.  The 
%       eval statements in aowfsGetFrameInfo below may also need 
%       modification in such circumstances.
%
% There are separately named variables for each frame: "xspotpos0000003", 
% etc. frmstrs{ii} contains the seven final numeric characters, which 
% specify the frame number, ii. Because frame 1 data is specified by 
% '0000000',  frame 224 is specified by '0000223', etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: aowfsGetFrameInfo.m 3051 2010-10-01 20:33:26Z amoran $

%% BEGIN_CODE

load(['slp',irigin]);
load(['int',irigin]);
numframes = nframes;
frmstrs = aowfsGetFrmstrs(numframes);
nsub = procInfo.nsubaps;
countstorad = procInfo.countstorad;
ivtocn2ryt = procInfo.ivtocn2ryt;  % parameter for cn2 calculation
apermask = procInfo.apermask;

for ii=1:1:numframes
   eval(['xpos=xsp',frmstrs{ii},';'])
   eval(['ypos=ysp',frmstrs{ii},';'])
   eval(['inten=int',frmstrs{ii},';'])
%   frameInfo(ii).time = timestamp;   % ENABLE THIS WHEN TIME INFO IS GIVEN
   slopex = zeros(nsub,1);
   slopey = zeros(nsub,1);
   irr = zeros(nsub,1);
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

%  compute irradiance quantities:

   irrBar = mean(irr);
   irrVar = (std(irr))^2;
   irrFNorm = irr/irrBar;
   irrFNVar = (std(irrFNorm))^2;
   cn2FN = ivtocn2ryt * irrFNVar;

   frameInfo(ii).irr = irr;
   frameInfo(ii).irrBar = irrBar;
   frameInfo(ii).irrVar = irrVar;
   frameInfo(ii).irrFNVar = irrFNVar;
   frameInfo(ii).cn2FN = cn2FN;

   if (nz ~= 0)
%  compute Zernikes

       slopes = [slopex' slopey'];
       zern = procInfo.winv*slopes.';
       frameInfo(ii).zern = zern;
   end


   if (flagR0slp ~= 0)
%  compute tilt structure functions

       [dxxx,dxxy,dyyx,dyyy,dxxxy,dyyxy]=mkslpstfb1(slopex,slopey,...
                                                    procInfo.xsub,procInfo.ysub,...
                                                    procInfo.xcenter,procInfo.ycenter,...
                                                    procInfo.spacing,procInfo.inradius);

       frameInfo(ii).dxxx = dxxx;
       frameInfo(ii).dxxy = dxxy;
       frameInfo(ii).dyyx = dyyx;
       frameInfo(ii).dyyy = dyyy;
       frameInfo(ii).dxxxy = dxxxy;
       frameInfo(ii).dyyxy = dyyxy;
   end
end

