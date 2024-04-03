function [opdiffull, opdifinv, xopd, yopd] = fullopdinv(opdif, opdnx, opdny, opddx, opddy, xact, yact)
% SYNTAX:
% [opdiffull, opdifinv, xopd, yopd] = fullopdinv(opdif, opdnx, opdny, ...
%                                     opddx, opddy, xact, yact)
% [opdiffull, opdifinv, xopd, yopd] = fullopdinv(opdif(imas,:), opdnx, ...
%                              opdny, opddx, opddy, xact(imas), yact(imas))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Builds the full OPD influence function matrix from the box-stored matrix. 
% Unnecessary OPD grid points are removed & grid point locations are 
% calculated in output space. The psuedo-inverse of full OPD influence 
% function matrix is then computed. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% opdif [ ] = The nact x nopd^2 box-stored OPD influence function matrix. 
%             nopd is size of the OPD grid influenced by one actuator and 
%             nact is the number of (master) actuators.
% opdnx, opdny [ ] = The size of the full OPD grid.
% opddx, opddy [ ] = Spacing of OPD grid.
% xact, yact [ ] = Actuator locations in output space.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% opdiffull [ ] = The full OPD influence matrix.
% opdifinv [ ] = The pseudo-inverse dimensioned the same as the transpose 
%                of opdiffull.
% xopd, yopd [ ] = Locations of necessary OPD grid points
% Notes:
%    The typical process for creating a reconstructor is:
%    1. Run AOGeom to specify the adaptive optics geometry.
%    2. Run AOInf to create the influence functions.
%    3. Load the resulting file into a clean Matlab workspace.
%    4. [opdiffull, opdifinv, xopd, yopd] = fullopdinv(opdif, opdnx, opdny, opddx, opddy, xact, yact)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
% AUTHOR: Keith Beardmore
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: fullopdinv.m 3468 2012-05-17 18:41:01Z keith $

%% BEGIN_CODE

msg = nargchk(7,7,nargin);
if (~isempty(msg))
   warning(['fullopdinv: ', msg]);
   help fullopdinv;
   return;
end;
%
if any(isnan(opdif))
   warning(['fullopdinv: ', 'opdif contains NaNs']);
   help fullopdinv;
   return;
end
[nact,nopdxy]=size(opdif);
disp(['OPD box-stored as nact x nopd^2 = ',int2str(nact),' x ',int2str(nopdxy)])
nopdxy=sqrt(nopdxy);
opdif=reshape(opdif,nact,nopdxy,nopdxy);
disp(['OPD re-shaped to nact x nopd x nopd = ',int2str(nact),' x ',int2str(nopdxy),' x ',int2str(nopdxy)])
%
disp(['Full OPD grid size is opdnx x opdny = ',int2str(opdnx),' x ',int2str(opdny)])
disp(['Full OPD grid spacing is opddx x opddy = ',num2str(opddx),' x ',num2str(opddy)])
%
actdx=min(diff(unique(xact)));
actnx=1+(max(xact)-min(xact))/actdx;
disp(['Actuator grid size = ',int2str(actnx)])
disp(['Actuator grid spacing = ',num2str(actdx)])
%
disp('Building Full OPD Influence Matrix ...')
xopd=nan*ones(opdnx,opdny);
yopd=nan*ones(opdnx,opdny);
opdiffull=nan*ones(nact,opdnx,opdny);
ihaf=floor(nopdxy/2);
rad2=max(xact.^2+yact.^2); % Maximum actuator radius
disp(['Maximum radius = ',num2str(sqrt(rad2)/opddx),' opd points or ',num2str(sqrt(rad2)),' m'])
obs2=min(xact.^2+yact.^2); % Minimum actuator radius
if(obs2>1e-10)
    disp(['Annular DM; Minimum radius = ',num2str(sqrt(obs2)/opddx),' opd points or ',num2str(sqrt(obs2)),' m'])
end
obs2=obs2-0.25*opddx^2-0.25*opddy^2; % Allow 1/2 OPD spacing for rounding error
rad2=rad2+0.25*opddx^2+0.25*opddy^2;
for i=1:nact
    % Actuator locations in OPD pixels
    % Actuator locations in OPD pixels
    ix=1+round(xact(i)/opddx)+floor(opdnx/2);
    iy=1+round(yact(i)/opddy)+floor(opdny/2);
    % OPD limits for this actuator in OPD pixels
    ilow=max(1,ix-ihaf);
    ihih=min(opdnx,ix+ihaf);
    jlow=max(1,iy-ihaf);
    jhih=min(opdny,iy+ihaf);
    for j=ilow:ihih
        % OPD pixel location in compressed matrix
        j1=j - ix + ihaf+1;
        for k=jlow:jhih
            % OPD pixel location in compressed matrix
            k1=k - iy + ihaf+1;
            % Check we're on the DM or within a subaperture
            x=(j-1-floor(opdnx/2))*opddx;
            y=(k-1-floor(opdny/2))*opddy;
            if(((x^2+y^2)<=rad2) && ((x^2+y^2)>=obs2))
                % Full OPD influence function matrix
                opdiffull(i,j,k)=opdif(i,j1,k1);
                if isnan(opdif(i,j1,k1))
                    disp(opdif(i,j1,k1))
                end
                % OPD pixel location in output space
                xopd(j,k)=x;
                yopd(j,k)=y;
            end
        end
    end
end
disp('Built Full OPD Influence Matrix')
%
% Re-shape OPD matrix
%disp(['Size(opdiffull)=',int2str(size(opdiffull))])
opdiffull=reshape(opdiffull,nact,opdnx*opdny);
%disp(['Size(opdiffull)=',int2str(size(opdiffull))])
% Remove untouched opd locations from OPD matrix
mask=(sum(isnan(opdiffull),1)~=nact);
opdiffull=opdiffull(:,mask);
mask=isnan(opdiffull);
opdiffull(mask)=0.0;
% Need to transpose
opdiffull=opdiffull';
disp(['Size of full OPD influence matrix is ',int2str(size(opdiffull))])
%
% Re-shape OPD locations
%disp(['Size(xopd)=',int2str(size(xopd))])
%disp(['Size(yopd)=',int2str(size(yopd))])
xopd=reshape(xopd,[],1);
yopd=reshape(yopd,[],1);
%disp(['Size(xopd)=',int2str(size(xopd))])
%disp(['Size(yopd)=',int2str(size(yopd))])
% Remove untouched opd locations
mask=(~isnan(xopd));
xopd=xopd(mask);
yopd=yopd(mask);
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
