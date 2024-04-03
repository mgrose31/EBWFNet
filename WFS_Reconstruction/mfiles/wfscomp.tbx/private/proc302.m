%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: proc302.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE
 

path('/nfs/u0/ablact/MODBC/src',path);
path('/nfs/u0/ablact/miscm',path);
path('/nfs/u0/ablact/aowfs',path);

outdir = [pwd,'/302am'];
datadir = '/nfs/u0/overflow/aowfs/incoming/302am'

cd(datadir);
amfiles = dir('slp*')
n = length(amfiles);

nz = 20;  % number of Zernike terms (nz=0 selects NONE)
flagR0z = 1; % Calculate r0 from Zernikes (flagR0z=0 selects NO CALCULATION)
flagR0slp = 1; % Calculate r0 from slope structure functions (flagR0slp=0 selects NO CALCULATION)
procInfo = aowfsGetProcInfo(nz);

for ii=1:n
   cd(datadir);
   irigin = amfiles(ii).name(4:14)
   ii
   tic;
   frameInfo = aowfsGetFrameInfo4(irigin,procInfo,nz,flagR0slp);
   'frameInfo', timing(1) = toc
   tic;
   accum = aowfsAccumFrame(frameInfo);
   'accum', timing(2) = toc
   tic;
   stats = aowfsComputeStats(accum,frameInfo,procInfo,nz,flagR0z,flagR0slp);
   'stats', timing(3) = toc
   tic;
   cd(outdir);
   save(['res',irigin], 'frameInfo','accum','stats','timing');
   'save', timing(4) = toc
   clear('frameInfo','accum','stats','timing');
end;

outdir = [pwd,'/302pm'];
datadir = '/nfs/u0/overflow/aowfs/incoming/302pm'

cd(datadir);
amfiles = dir('slp*')
n = length(amfiles);

nz = 20;  % number of Zernike terms (nz=0 selects NONE)
flagR0z = 1; % Calculate r0 from Zernikes (flagR0z=0 selects NO CALCULATION)
flagR0slp = 1; % Calculate r0 from slope structure functions (flagR0slp=0 selects NO CALCULATION)
procInfo = aowfsGetProcInfo(nz);

for ii=1:n
   cd(datadir);
   irigin = amfiles(ii).name(4:14)
   ii
   tic;
   frameInfo = aowfsGetFrameInfo4(irigin,procInfo,nz,flagR0slp);
   'frameInfo', timing(1) = toc
   tic;
   accum = aowfsAccumFrame(frameInfo,nz,flagR0slp);
   'accum', timing(2) = toc
   tic;
   stats = aowfsComputeStats(accum,frameInfo,procInfo,nz,flagR0z,flagR0slp);
   'stats', timing(3) = toc
   tic;
   cd(outdir);
   save(['res',irigin], 'frameInfo','accum','stats','timing');
   'save', timing(4) = toc
   clear('frameInfo','accum','stats','timing');
end;

