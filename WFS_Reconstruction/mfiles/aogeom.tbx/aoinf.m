function aoinf(infile, outfile, nxopd, dxopd)
% SYNTAX:
% aoinf(infile, outfile, nxopd, dxopd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% aoinf calls the external program AOInf in order to compute the 
% deformable mirror OPD and slope influence functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% infile [ ] = is the name of input file from which the AO geometry 
%              specifications are to be read. Usually, but not necessarily,
%              this file is created by AOGeom. The filename specification 
%              must include the .mat suffix.
% outfile [ ] = is the name of file to which the resulting combined
%               reconstructors and AO geometry are to be written. This 
%               file is then loaded into Matlab for computation of the 
%               reconstructor using AORecon or specified as a constructor 
%               argument to AOModel. The filename specification must 
%               include the .mat suffix.
% nxopd [ ] = Specifies the number of points across the square grid
%             superimposed on the DM for the purpose of resolving surface 
%             deflections on the OPD mesh. This number does not need to be 
%             as large as the propagation grid because the OPD influence 
%             function needs only to provide enough points to cover the 
%             entire mirror. The OPD mesh that is used is one at which the 
%             origin lies precisely at nxopd/2 + 1 point. Hence, the mesh 
%             is not symmetric and one must be careful to use nxopd and 
%             dxopd combination which covers the entire physical region of 
%             interest.
% dxopd [ ] = Specifies the physical spacing between grid points of the 
%             square grid superimposed on the DM for the purpose of 
%             resolving surface deflections on the OPD mesh. This 
%             parameter needs to be small enough to properly model the 
%             shape of the mirror and must to be smaller than  the spacing 
%             between actuators. It is usually sufficient to use six to 
%             twelve points per actuator spacing. Also, the actuator 
%             spacing must  be an even multiple of dxopd.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% None. A series of messages are printed which indicate the success or
% lack thereof in computing the influence functions and writing the file.
% NOTES:
% The two main reasons that AOInf may not complete successfully is because
% the user has specified a dxopd which does not go into the actuator 
% spacing an integral number of times. Or the nxopd, dxopd combination 
% does not entirely cover the grid of actuators to be covered. The best 
% way to avoid these problems is to specify dxopd as an integer divisor of 
% the actuator spacing (i.e., dxopd = dxact/n, where n is an integer). And 
% use the command opdc to generate opd coordinates and plot up with xact, 
% yact, xsub, and ysub which appear in the AOGeom output files.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: aoinf.m 4212 2023-05-05 17:52:03Z jtellez $

%% BEGIN_CODE

msg = nargchk(4,4,nargin);
if (~isempty(msg))
   warning(['aoinf: ', msg]);
   help aoinf;
   return;
end;
if (~ischar(infile))
   warning(['aoinf: infile must be a character string containing the input file name.']);
   help aoinf;
   return;
end;
if (~ischar(outfile))
   warning(['aoinf: outfile must be a character string containing the input file name.']);
   help aoinf;
   return;
end;
if (~isempty(nxopd) & isreal(nxopd))
   if ((nxopd <= 0) & (nxopd ~= fix(nxopd)))
      warning(['aoinf: nxopd must be a positive integer.']);
      help aoinf;
      return;
   end;
else
   warning(['aoinf: nxopd must be a positive integer.']);
   help aoinf;
   return;
end;
if (~isempty(dxopd) & isreal(dxopd))
   if (dxopd <= 0)
      warning(['aoinf: dxopd must be positive.']);
      help aoinf;
      return;
   end;
else
   warning(['aoinf: dxopd must be positive.']);
   help aoinf;
   return;
end;
wtdir = getenv('WT_DIR');
cmdstr = ['!"' wtdir '"\bin\AOInf ', infile, ' ', outfile, ' ', num2str(fix(nxopd)), ' ', num2str(dxopd,'%20.14g')];
eval(cmdstr);