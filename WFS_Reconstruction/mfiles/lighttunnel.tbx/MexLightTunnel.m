function [img,tmxh] = MexLightTunnel(imgp, p, tx, ty, bx, by, rho, tmxh, keepmex, FL, z)
% SYNTAX:
% [img] = MexLightTunnel(imgp, p, tx, ty, bx, by, rho, tmxh, keepmex, FL, z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Adds specified tilt and blur effects to input image.
% Uses mex system of WaveTrain LightTunnel code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% imgp.g [matrix] = Input image at reflector plane (must be double)
% p/p.g [scalar/matrix] = power multiplier of output image values
% tx.g [matrix] = X-translations for image pixel locations (m)
% ty.g [matrix] = Y-translations for image pixel locations (m)
% bx/bx.g [scalar/matrix] = Row blur function per pixel (m, default 0)
% by/by.g [scalar/matrix] = Column blur function per pixel (m, default 0)
%    Note: Tilts and blurs are in reflector plane units,
%    e.g., tx is apparent translation of PS due to tilt as seen from camera
%    rho/rho.g [scalar/matrix] = PSF orientation per pixel (-1 to 1, def 0)
% For matrix inputs:
% *.x [vector] = X position of image pixels, increasing down columns (arb)
% *.y [vector] = Y position of image pixels, increasing across rows (arb)
%    Note: All inputs grids must be defined in the reflector plane
%    and be the same size and spacing
% tmxh [handle] = handle of a previously created mex system instance to be 
%                 used, default is [] so a new instance is created.
% keepmex [ ] = keep current mex system instance (true/false, default is false)
% FL [ ] = camera focal length to calculate magnification (m, default = 1.0)
% z [ ] = propagation distance to calculate magnification (m, default = 1.0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% img.g [matrix] = imgp with applicable tilt and blur added image
%                (camera plane image)
% img.x [vector] = Image x-coordinates (increasing down row)
% img.y [vector] = Image y-coordinates (increasing across column)
%              (camera pixel coordinates)
% tmxh [handle] = handle of current mex system instance, [] if none
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mex usage note:
% The matlab path must include directories for both this .m file and for
% the mex .dll (or .cpp file to build the .dll)
% Simple usage: set tmxh = [] & keepmex = false (default values)
%               A new mex instance is created & destroyed each time
% Efficient usage: first: tmxh = [] & keepmex = true
%                  other: tmxh = previous tmxh output & keepmex = true
%                  final: tmxh = previous tmxh output & keepmex = false
%                         The same mex instance is used each time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Keith m Beardmore
% (c) 2007 MZA Associates Corporation, Albuquerque, NM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: MexLightTunnel.m 3036 2010-09-23 21:22:50Z amoran $

%% BEGIN_CODE

% Check that minimum inputs are defined
if nargin < 4
    error('At least four inputs are required; imgp, p, tx, & ty')
else
    if nargin < 5
        bx = [];
    end
    if nargin < 6
        by = [];
    end
    if nargin < 7
        rho = [];
    end
end

% Turn scalars into grids
if ~isstruct(p)
    tmp=p;p=[];
    p.g=repmat(tmp,size(imgp.g));
    clear tmp
    p.x=imgp.x;
    p.y=imgp.y;
end
if isempty(bx), bx = 0.0; end
if ~isstruct(bx)
    tmp=bx;bx=[];
    bx.g=repmat(tmp,size(imgp.g));
    clear tmp
    bx.x=imgp.x;
    bx.y=imgp.y;
end
if isempty(by), by = 0.0; end
if ~isstruct(by)
    tmp=by;by=[];
    by.g=repmat(tmp,size(imgp.g));
    clear tmp
    by.x=imgp.x;
    by.y=imgp.y;
end
if isempty(rho), rho = 0.0; end
if ~isstruct(rho)
    tmp=rho;rho=[];
    rho.g=repmat(tmp,size(imgp.g));
    clear tmp
    rho.x=imgp.x;
    rho.y=imgp.y;
end

% Check that all grids have the same dimensions
if ~(isequal(imgp.x, p.x, tx.x, ty.x, bx.x, by.x, rho.x) ...
  && isequal(imgp.y, p.y, tx.y, ty.y, bx.y, by.y, rho.y))
    error('All input grids must have same dimensions')
end

mexname='LightTunnelingRunMex';
dllname=['mex',mexname,'.dll'];
cppname=['mex',mexname,'.cpp'];

% Make sure mex system is built
if ~exist(dllname,'file')
    if exist(cppname,'file')
        rundir=pwd;
        mexdir=which(cppname);mexdir=mexdir(1:end-length(cppname));
        disp(['Building ',dllname,' in directory ',mexdir,' (moving from ',rundir,')'])
        cd(mexdir)
        tmxmkdll(mexname)
        cd(rundir)
        disp(['Back in directory ',rundir])
    else
        error(['Cannot find ',cppname,' to build ',dllname])
    end
end

% Create instance of dll system if it does not already exist
if ~exist('tmxh','var'), tmxh = []; end
if isempty(tmxh)
    disp(['Creating ',mexname,' instance'])
    % Change mex system parameters
    if ~exist('FL','var'), FL=[]; end
    if isempty(FL), FL = 1.0; end
    if ~exist('z','var'), z=[]; end
    if isempty(z), z = 1.0; end
    params = tmxparams(mexname);
    params.FL = FL;
    params.z = z;
    tmxh = tmxcreate(mexname, params);
    tmxh = tmxrecoff(tmxh); % Turn off all .trf recording
end

% Pass input values
tmxh = tmxsetinp(tmxh, 'reflectedIntensity', imgp); % Returned intensity at reflector plane
tmxh = tmxsetinp(tmxh, 'relativeIntensity', p); % Power in PSF relative to PS power
tmxh = tmxsetinp(tmxh, 'xDisplacement', tx); % PSF translation in reflector plane units
tmxh = tmxsetinp(tmxh, 'yDisplacement', ty); % PSF translation in reflector plane units
tmxh = tmxsetinp(tmxh, 'xSpread', bx); % PSF standard deviation in reflector plane units
tmxh = tmxsetinp(tmxh, 'ySpread', by); % PSF standard deviation in reflector plane units
tmxh = tmxsetinp(tmxh, 'orientation', rho); % PSF orientation

% Run mex system
tmxh = tmxrun(tmxh, tmxh.stoptime+0.001);

% Get Light-Tunneled output
img = tmxload(tmxh, 'lighttunnel.fpaImage');

% Delete the mex instance if requested
if ~exist('keepmex','var'), keepmex=[];end
if isempty(keepmex), keepmex = false; end
if ~keepmex
    disp(['Closing ',mexname,' instance'])
    trfname=tmxh.trfname;
    tmxh = tmxdel(tmxh);
    if exist(trfname,'file'), delete(trfname); end
end
