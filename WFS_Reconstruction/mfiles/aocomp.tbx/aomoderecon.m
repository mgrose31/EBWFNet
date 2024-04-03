function [recon, rsv, ru, rs, rv]  = aomodrecon(mvtos, nsingin)
% SYNTAX: 
% recon = aomoderecon(mvtos, nsingin)
% rsv = aomoderecon(mvtos)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Compute a least-squares reconstructor from a mode constrained slope 
% influence matrix. The reconstructor can have a user-specified number of 
% singular values removed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% mvtos [ ] = The 2*nsub x nmod mode constrained slope influence function 
%             matrix. nsub is the  number of subapertures and nmod is the 
%             number of modes. This is the usual slope influence matrix, 
%             post multiplied by a matrix that generates actuator 
%             displacments from mode coefficients. The rows must be 
%             arranged such that nsub x-slopes are followed by nsub 
%             y-slopes.
% nsing [ ] = The number of singular values to be removed. When not 
%             specified this routine does not compute the reconstructor. 
%             Rather it returns a vector containing the singular values 
%             for the purpose of deciding how many to remove on a 
%             subsequent call to this routine.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% recon [ ] = The reconstructor dimensioned the same as the transpose of 
%             mvtos.
% rsv [ ] = A vector containing the singular values of mvtos.
% Notes:
%    The typical process for creating a reconstructor is:
%    1. Run AOGeom to specify the adaptive optics geometry.
%    2. Run AOInf to create the influence functions.
%    3. Load the resulting file into a clean Matlab workspace.
%    4a. Generate the mode-actuator matrix & post multiply mvtos
%    4. Execute rsv=aomoderecon(mvtos)
%    5. Inspect rsv. Usually there is a clear breakpoint where the lowest
%       few singular values are much less than all the others. It is often 
%       good to project these out. Determine the number to be suppressed 
%       and set the variable ns to it.
%    6. Execute recon=aomoderecon(mvtos,ns) to create the reconstructor.
%    7. Pre multiply recon by the mode-actuator matrix
%    8. Save the workspace. The file can be used as input to WaveTrain 
%       (usually as an argument to a TasatDMModel constructor).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: aomoderecon.m 3027 2010-09-21 21:04:10Z amoran $

%% BEGIN_CODE

msg = nargchk(1,2,nargin);
if (~isempty(msg))
   warning(['aomoderecon: ', msg]);
   help aomoderecon;
   return;
end;
nsub2 = size(mvtos, 1);
if (2*floor(nsub2/2) ~= nsub2)
   warning(['aomoderecon: The inner dimension of mvtos must be even.']);
   help aomoderecon;
   return;
else
   nsub = nsub2/2;
end;
%
nsing = [];
if (nargin > 1)
   if (~isempty(nsingin))
      nsing = nsingin;
   end;
end;
%
[ru,rs,rv]=svd(mvtos);
%
rsv=diag(rs);
rmaxsv=max(rsv);
rnsv=size(rsv,1);
rsv1=ones(rnsv,1)./rsv;
if (isempty(nsing))
   recon=rsv;
   fprintf(1,'The modal reconstructor WAS NOT computed.\n');
   return;
end;
if (nsing > 0)
   rsv1((rnsv-nsing+1):rnsv)=zeros(nsing,1);
end;
rsn1=[diag(rsv1),zeros(size(rsv1,1),size(ru,1)-size(rsv1,1))];
recon=rv*rsn1*ru';
fprintf(1,'The modal reconstructor WAS computed.\n');

