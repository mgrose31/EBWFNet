function [opdifinv, sv] = opdinv(opdif, nsingin)
% SYNTAX:
% [opdifinv, sv] = opdinv(opdif, nsingin)
% sv = opdinv(opdif)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Compute the psuedo-inverse of the OPD influence function matrix. The 
% inverse can have a user-specified number of singular values removed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% opdif [ ] = The nopd x nact OPD influence function matrix. nopd is the 
%             number of OPD sample points and nact is the number of 
%             (master) actuators.
% nsingin [ ] = The number of singular values to be removed. When not 
%               specified this routine does not compute the inverse. 
%               Rather it returns a vector containing the singular values 
%               for the purpose of deciding how many to remove on a 
%               subsequent call to this routine.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% opdifinv [ ] = The pseudo-inverse dimensioned the same as the transpose 
%                of opdif.
% sv [ ] = A vector containing the singular values of opdif.
% Notes:
%    The typical process for creating opdifinv is:
%    1. Run AOGeom to specify the adaptive optics geometry.
%    2. Run AOInf to create the influence functions.
%    3. Load the resulting file into a clean Matlab workspace.
%    4. Execute sv=opdinv(opdif)
%    5. Inspect sv. There probably won't be any singular modes to remove.
%       Usually there is a clear breakpoint where the lowest
%       few singular values are much less than all the others. It is often 
%       good to project these out. Determine the number to be suppressed 
%       and set the variable ns to it.
%    6. Execute opdifinv = opdinv(opdif, nsingin) to create the pseudo-inverse.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
% AUTHOR:Keith@MZA 3/22/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: opdinv.m 3027 2010-09-21 21:04:10Z amoran $

%% BEGIN_CODE

msg = nargchk(1,2,nargin);
if (~isempty(msg))
   warning(['opdinv: ', msg]);
   help opdinv;
   return;
end;
%
nsing = [];
if (nargin > 1)
   if (~isempty(nsingin))
      nsing = nsingin;
   end;
end;
%
[u,s,v]=svd(opdif);
sv=diag(s);
clear s
maxsv=max(sv);
nsv=size(sv,1);
sv1=ones(nsv,1)./sv;
%
if (isempty(nsing))
   opdifinv=sv;
   fprintf(1,'The inverse OPD influence function matrix WAS NOT computed.\n');
   return;
end;
if (nsing > 0)
   sv1((nsv-nsing+1):nsv)=zeros(nsing,1);
end;
sn1=[diag(sv1),zeros(size(sv1,1),size(u,1)-size(sv1,1))];
opdifinv=v*sn1*u';
fprintf(1,'The inverse OPD influence function matrix WAS computed.\n');
