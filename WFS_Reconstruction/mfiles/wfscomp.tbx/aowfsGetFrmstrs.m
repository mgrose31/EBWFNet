function frmstrs = aowfsGetFrmstrs(numframes);
% SYNTAX: 
% frmstrs = aowfsGetFrmstrs(numframes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:  
% aowfsGetFrmstrs produces a six-character string which is required to 
% access aowfs data stored in variables named by frame number.  For 
% example, the  x spot position data for frame 224 is in variable 
% 'xspotpos0000223'.  The function  generates the strings required to read 
% each frame from a data set which contains the number of frames input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% numframes:  The total number of frames in the data set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% frmstrs [ ] =  A cell array of length numframes, eachcell of which 
%                contains the seven numeric final characters of the data 
%                variables for that frame. frmstrs{ii} is used in an eval 
%                statement by the function aowfsGetFrameInfo.
%
% NOTE: THIS FUNCTION MUST BE MODIFIED if the naming convention for the 
%       data variables is changed, and the eval statements in 
%       aowfsGetFrameInfo may need to be modified.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: aowfsGetFrmstrs.m 3051 2010-10-01 20:33:26Z amoran $

%% BEGIN_CODE
 
numstring = num2str(numframes);
ndigits = length(numstring);
nzeros = 7 - ndigits;
if (nzeros < 1)
  error('more than 1e6 frames')
end
zerstr{1} = '0';
zerstr{2} = '00';
zerstr{3} = '000';
zerstr{4} = '0000';
zerstr{5} = '00000';
zerstr{6} = '000000';
zerstr{7} = '0000000';
for ii=1:1:numframes
  jj = ii-1;
  digstr = num2str(jj);
  lndig = length(digstr);
  lnzer = 7 - lndig;
  frmstrs{ii}=[zerstr{lnzer},digstr];
end

