function [r0XX,r0XY,r0XXY,r0YX,r0YY,r0YXY,dxxx,dxxy,dyyx,dyyy,dxxxy,dyyxy,dxxxth,dxxyth,dxxxyth]=r0slpstrf(s,scal);
% SYNTAX:
% [r0XX,r0XY,r0XXY,r0YX,r0YY,r0YXY,dxxx,dxxy,dyyx,dyyy,dxxxy,dyyxy, ...
%                                  dxxxth,dxxyth,dxxxyth]=r0slpstrf(s,scal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% computes the experimental & theoretical slope structure functions and
% compares the two in order to extract estimates of r0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% s [ ] = the structure containing slope and subaperture information
%         concerning the WFS frames that are to be analyzed
% scal [ ] = the structure containing 
% spar [ ] = the structure containing the system specific information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS: 
% r0XX , r0XY, r0XXY the estimates of r0 using X-slope information in
%        the X Y and XY(diagonal) directions
%r0YX , r0YY, r0YXY the estimates of r0 using Y-slope information in
%       the X Y and XY(diagonal) directions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: r0slpstrf.m 3051 2010-10-01 20:33:26Z amoran $

%% BEGIN_CODE
 
ndatasets = size(s.xslp,2);

d = scal.spacing;

option=foptions;

xsub=scal(1).xsub;
ysub=scal(1).ysub;

[dxxx,dxxy,dyyx,dyyy,dxxxy,dyyxy]=mkslpstfb1(s.xslp,s.yslp,xsub,ysub,scal.xcenter,scal.ycenter,scal.spacing,scal.inradius);

if (max(size(dxxx)) == 1)
   r0XX=0;
   r0XY=0;
   r0XXY=0;
   r0YX=0;
   r0YY=0;
   r0YXY=0;
   dxxxth=0;
   dxxyth=0; 
   dxxxyth=0;
   return
end;

npts=length(dxxx);
npts1=length(dxxy);
if ((npts > 3) & (npts1 > 3))
   [dxxxth,dxxyth,dxxxyth]=cpaperstrt(npts);
   tmp1=[1:1:fix(2*npts/3)];
   tmp1a=[1:1:fix(2*npts1/3)];
else
   'Insufficient Number of Points'
end

bestx=fmin('trialfun',.01,6000,option,dxxx(tmp1),dxxxth(tmp1));
r0XX = (bestx*(d^(1/3)))^(-3/5);
besty=fmin('trialfun',.01,6000,option,dxxy(tmp1a),dxxyth(tmp1a));
r0XY = (besty*(d^(1/3)))^(-3/5);

bestx=fmin('trialfun',.01,6000,option,dyyx(tmp1),dxxyth(tmp1));
r0YX = (bestx*(d^(1/3)))^(-3/5);
besty=fmin('trialfun',.01,6000,option,dyyy(tmp1a),dxxxth(tmp1a));
r0YY = (besty*(d^(1/3)))^(-3/5);

npts=length(dxxxy);
npts1=length(dyyxy);
if npts>3 & npts1>3

   [dxxxth,dxxyth,dxxxyth]=cpaperstrt(npts);
   tmp1=[1:1:fix(2*npts/3)];
   tmp1a=[1:1:fix(2*npts1/3)];
else
   'Insufficient Number of Points'
end

bestx=fmin('trialfun',.01,6000,option,dxxxy(tmp1),dxxxyth(tmp1));
r0XXY = (bestx*(d^(1/3)))^(-3/5);
besty=fmin('trialfun',.01,6000,option,dyyxy(tmp1a),dxxxyth(tmp1a));
r0YXY = (besty*(d^(1/3)))^(-3/5);

