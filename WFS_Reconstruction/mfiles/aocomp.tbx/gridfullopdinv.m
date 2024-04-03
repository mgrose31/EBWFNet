function [opdiffull, opdifinv, xopd, yopd] = gridfullopdinv(actinfl, xact, yact)
% SYNTAX:
% [opdiffull, opdifinv, xopd, yopd] = gridfullopdinv(actinfl, xact, yact)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Builds the full OPD influence function matrix from the
% A set of influence function grids.
% Unnecessary OPD grid points are removed & grid point
% locations are calculated in output space. The psuedo-inverse of full OPD
% influence function matrix is then computed. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% actinfl [ ] = Cell array where actinfl{i}.x, actinfl{i}.y, actinfl{i}.g 
%               is the influence function for actuator i
% xact, yact [ ] = Actuator locations in output space.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% opdiffull [ ] = The full OPD influence matrix.
% opdifinv [ ] = The pseudo-inverse dimensioned the same as the transpose 
%                of opdiffull.
% xopd, yopd [ ] = Locations of necessary OPD grid points
%
% Notes:
% 1. Assumes all influence functions are the same size, span entire DM 
% and do not need to be translated
% 2. The output OPD grid points do not necessarily relate to the actuator
% influence function mesh:
%    (a) The OPD grid may be coarser (downsampled influence functions) in
%        order to reduce memory requirements
%    (b) The OPD points may be recalculated (interpolated influence functions)
%        in order to force alignment with actuators and suabapertures.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
% AUTHOR: Keith Beardmore
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: gridfullopdinv.m 3141 2010-11-11 21:06:28Z keith $

%% BEGIN_CODE

msg = nargchk(3,3,nargin);
if (~isempty(msg))
   warning(['gridfullopdinv: ', msg]);
   help gridfullopdinv;
   return;
end;
%
nact=length(actinfl);
disp(['OPD stored as nact [= ',int2str(nact),'] grids '])
%
opdx = actinfl{1}.x;
opdy = actinfl{1}.y;
opdnx = length(opdx);
opdny = length(opdy);
opddx = mean(diff(opdx));
opddy = mean(diff(opdy));
disp(['Full OPD grid size is opdnx x opdny = ',int2str(opdnx),' x ',int2str(opdny)])
disp(['Full OPD grid spacing is opddx x opddy = ',num2str(opddx),' x ',num2str(opddy)])
%
% Test that all actuator locations coincide with OPD mesh
dactopd2 = 0;
interp=false;
for i=1:nact
    dactopd2 = max(dactopd2, min((xact(i)-actinfl{1}.x).^2 + (yact(i)-actinfl{1}.y).^2));
end
if sqrt(dactopd2) > min(opddx,opddy)/1000;
    disp('gridfullopdinv.m: OPD grid does not subsample actuator locations, interpolating OPD meshes')
    interp=true;
else % Test to  see if we need to downsample
    % Will probably have to downsample to get a small enough OPD matrix to invert
    nact1D = 2*sqrt(nact/pi); % Estimate number of actuators across
    nopdxy = 6*nact1D; % Approximate number of OPD points across (6 per actuator)
    if round(max(opdnx,opdny)/nopdxy) > 1
        disp('gridfullopdinv.m: OPD grid is too large, downsampling OPD meshes')
        interp=true;
    end
end
clear dactopd2
%
if interp
    opddx = mean(diff(unique(xact)))/5;% 6 OPD points per actuator
    opddy = mean(diff(unique(yact)))/5;
    opdx = min(xact):opddx:max(xact);
    opdy = min(yact):opddy:max(yact);
    opdnx = length(opdx);
    opdny = length(opdy);
    disp(['Revised OPD grid size is opdnx x opdny = ',int2str(opdnx),' x ',int2str(opdny)])
    [xx,yy]=meshgrid(actinfl{1}.x,actinfl{1}.y);
    for i=1:nact
        actinfl{i}.g=interp2(xx,yy,actinfl{i}.g,opdx,opdy','spline');
    end
    clear xx yy
end
clear interp
%
disp('Building Full OPD Influence Matrix ...')
rad2=max(xact.^2+yact.^2); % Maximum actuator radius
disp(['Maximum radius = ',num2str(sqrt(rad2)/opddx),' opd points or ',num2str(sqrt(rad2)),' m'])
obs2=min(xact.^2+yact.^2); % Minimum actuator radius
if(obs2>1e-10)
    disp(['Annular DM; Minimum radius = ',num2str(sqrt(obs2)/opddx),' opd points or ',num2str(sqrt(obs2)),' m'])
end
obs2=obs2-0.25*opddx^2-0.25*opddy^2; % Allow 1/2 OPD spacing for rounding error
rad2=rad2+0.25*opddx^2+0.25*opddy^2;
opdiffull=nan*ones(nact,opdnx,opdny);
[xopd,yopd] = meshgrid(opdx,opdy);
r2 = xopd.^2+yopd.^2;
mask = r2<rad2 & r2>obs2;
for i=1:nact
    opdiffull(i,:,:)=mask.*actinfl{i}.g';
end

disp('Built Full OPD Influence Matrix')
%
% Re-shape OPD matrix & OPD locations
opdiffull=reshape(opdiffull,nact,opdnx*opdny);
xopd=reshape(xopd,[],1);
yopd=reshape(yopd,[],1);
% Remove untouched opd locations from OPD matrix
mask=(sum((opdiffull==0),1)~=nact);
opdiffull=opdiffull(:,mask);
xopd=xopd(mask);
yopd=yopd(mask);
% Need to transpose
opdiffull=opdiffull';
disp(['Size of full OPD influence matrix is ',int2str(size(opdiffull))])
disp(['Size of X OPD coordinate vector is ',int2str(size(xopd))])
disp(['Size of Y OPD coordinate vector is ',int2str(size(yopd))])
%
disp('Inverting Full OPD Influence Matrix ...')
sv=opdinv(opdiffull);
figure
semilogy(sv,'o-')
title('Singular values for full OPD influence function matrix')
%ns=input('How many singular modes?');
ns=0;disp('Assuming no singular modes')
opdifinv=opdinv(opdiffull,ns);    % Compute the opd pseudo-inverse.
disp('Inverted Full OPD Influence Matrix')
