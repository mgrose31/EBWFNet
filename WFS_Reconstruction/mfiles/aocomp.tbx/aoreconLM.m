function [recon, rsv] = aoreconLM(mvtos, nsingin, remtiltin, nsubin)
% 'aorecon' modified to minimize memory overhead in order to build
% reconstructors for LARGE (2000+ actuator) systems - Keith - March 2005
% Note: to minimize memory overhead, load only 'mvtos' from the .mat file
% into the workspace prior to calling this function.
% 1) ru,rs,rv no longer avialable as output - clear'd as soon as possible
% 2) rb no longer avialable as output - clear'd as soon as possible
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYNTAX:
% [recon, rsv] = aorecon(mvtos, nsing, {remtilt},  {nsubin})
% rsv = aorecon(mvtos)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Compute a least-squares reconstructor from a slope influence matrix. The 
% reconstructor can include tilt or project tilt out and can have a 
% user-specified number of singular values removed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% mvtos [ ] = The 2*nsub x nact slope influence function matrix. nsub is 
%             the number of subapertures and nact is the number of (master)
%             actuators. The rows must be arranged such that nsub x-slopes 
%             are followed by nsub y-slopes.
% nsing [ ] = The number of singular values to be removed. When not 
%             specified this routine does not compute the reconstructor. 
%             Rather it returns a vector containing the singular values 
%             for the purpose of deciding how many to remove on a 
%             subsequent call to this routine.
% remtilt [ ] = A flag which when non-zero indicates that tilt is to be
%               projected out. If not specified, tilt is projected out. To 
%               obtain a tilt-included reconstructor, specify zero for 
%               this argument.
% nsubin [ ] = Number of subapertures. Only necessary if it is not twice 
%              the first dimension of mvtos, e.g., if building the waffle 
%              constrained reconstructor.
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
%    4. Execute rsv=aorecon(mvtos)
%    5. Inspect rsv. Usually there is a clear breakpoint where the lowest
%       few singular values are much less than all the others. It is often 
%       good to project these out. Determine the number to be suppressed 
%       and set the variable ns to it.
%    6. Execute recon=aorecon(mvtos,ns) to create the tilt-removed 
%       reconstructor.
%    7. Execute reconti=aorecon(mvtos,ns,0) to create the tilt-included
%       reconstructor.
%    8. Save the workspace. The file can be used as input to WaveTrain 
%       (usually as an argument to a TasatDMModel constructor).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: aoreconLM.m 3027 2010-09-21 21:04:10Z amoran $

%% BEGIN_CODE

msg = nargchk(1,4,nargin);
if (~isempty(msg))
   warning(['aoreconLM: ', msg]);
   help aoreconLM;
   return;
end;
clear msg;pack
nsub2 = size(mvtos, 1);
if (nargin > 3)
   if (~isempty(nsubin))
      nsub2 = 2*nsubin;
   end;
   clear nsubin
end;
if (2*floor(nsub2/2) ~= nsub2)
   warning(['aoreconLM: The inner dimension of mvtos must be even.']);
   help aoreconLM;
   return;
else
   nsub = nsub2/2;
end;
clear nsub2;pack
remtilt = 1;
if (nargin > 2)
   if (~isempty(remtiltin))
      remtilt = remtiltin;
      clear remtiltin;pack
   end;
end;
nsing = [];
if (nargin > 1)
   if (~isempty(nsingin))
      nsing = nsingin;
      clear nsingin;pack
   end;
end;
%
disp('mvtos, nsing, nsub, remtilt')
whos
[ru,rs,rv]=svd(mvtos);
disp('mvtos, nsing, nsub, remtilt, rs, ru, rv')
whos
%
rsv=diag(rs);
clear rs;pack
if (isempty(nsing))
   recon=rsv;
   fprintf(1,'The reconstructor WAS NOT computed.\n');
   return;
end;
%
rnsv=size(rsv,1);
rsv1=ones(rnsv,1)./rsv;
clear rsv;pack
if (nsing > 0)
   rsv1((rnsv-nsing+1):rnsv)=zeros(nsing,1);
end;
clear rnsv nsing;pack
% Store 'mvtos' to disk
save tmp.mat mvtos
clear mvtos;pack
disp('nsub, remtilt, rsv1, ru, rv')
whos
%
rsn1=[diag(rsv1),zeros(size(rsv1,1),size(ru,1)-size(rsv1,1))];
clear rsv1;pack
recon=rv*rsn1*ru';
clear rv rsn1;pack
%
disp('nsub, recon, remtilt, ru')
whos
%
if (remtilt)
   rsx=[ones(nsub,1);zeros(nsub,1)];
   rsy=[zeros(nsub,1);ones(nsub,1)];
   if (size(ru,1) > (2*nsub))
      rsx = [rsx;zeros(nsub,1)];
      rsy = [rsy;zeros(nsub,1)];
   end;
   clear ru;pack
% Store 'recon', 'rsx', 'rsy' to disk
   save tmp.mat recon rsx rsy -append
   clear recon rsx rsy;pack
% Load 'mvtos' from disk
   load tmp.mat mvtos
   disp('mvtos, nsub, remtilt')
   whos
   rb=pinvLM(mvtos);
   clear mvtos;pack
% Load 'rsx' & 'rsy' from disk
   load tmp.mat rsx rsy
   rvx=rb*rsx;
   rvy=rb*rsy;
   clear rb;pack
   rq=orth([rvx,rvy]);
   clear rvx rvy;pack
   rp=eye(size(rq,1))-(rq*rq');
   clear rq;pack
   % Store 'rp' to disk
   save tmp.mat rp -append
   clear rp;pack
% Following line expanded to reduce memory requirements
   %rp1=eye(size(rsx,1))-((1/nsub)*(rsx*rsx' + rsy*rsy'));
   rp1=rsy*rsy';
   clear rsy;pack
   srsx=size(rsx,1);
   disp('nsub, rp1, rsx, remtilt')
   whos
   rp1=rp1+rsx*rsx';
   clear rsx;pack
   rp1=-(1/nsub)*rp1;
   clear nsub;pack
   rp1=rp1+eye(srsx);
   clear srsx;pack
   %
% Load 'recon' & 'rp' from disk
   load tmp.mat recon rp
   recon=rp*recon*rp1';
   disp('recon, remtilt, rp, rp1')
   whos
   clear rp rp1;pack
   fprintf(1,'The tilt-removed reconstructor WAS computed.\n');
else
   clear nsub ru;pack
   disp('recon, remtilt')
   whos
   fprintf(1,'The tilt-included reconstructor WAS computed.\n');
end;
%
clear remtilt;pack
disp('recon')
whos
