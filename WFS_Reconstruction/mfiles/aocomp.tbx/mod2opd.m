function [g] = mod2opd(xopd, yopd, orders_modes, rad, arad, ordering)
% SYNTAX:
% g = mod2opd(xopd, yopd, orders, rad)
% g = mod2opd(xopd, yopd, orders, rad, arad)
% g = mod2opd(xopd, yopd, modes, rad, [], ordering)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Compute the relationship between Zernike modes of the specified orders 
% and mirror displacement to give the correct OPD at locations across the DM.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% xopd [ ] = nlocs x 1 vector of x-coords of OPD locations
% yopd [ ] = nlocs x 1 vector of y-coords of OPD locations
% orders [ ] = The norder x 1 vector of orders of Zernike modes to be used,
%              e.g., 1= x & y tilt, 2 = focus, 0/90, & +/-45 astigmatism.
%              Note: nmod = Sum_{i=1:norder}(1+order(i))
% modes [] = The nmode x 1 vector of Zernike modes to be used. Requires the
%              'ordering' parameter to be set
% rad [ ] = Radius of the DM / beam/ aperture
% arad [ ] = Annular Radius (optional). If specified, annular Zernike 
%            polynomials are used. 
% ordering = Zernike ordering scheme (optional); 1=Noll, 2=Malacara, 3=Wyant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% g [ ] = The nloc x nmod matrix relating zernike coefficients to OPD. The
%         output is the mirror displacment, i.e 1/2 of the OPD
% Notes:
%    The typical process for creating a reconstructor is:
%    1. Run AOGeom to specify the adaptive optics geometry.
%    2. Run AOInf to create the influence functions.
%    3. Load the resulting file into a clean Matlab workspace.
%    4a. Generate the mode-opd & opd-actuator matrices & post multiply 
%        mvtos
%    4. Execute rsv=aorecon(mvtos)
%    5. Inspect rsv. Usually there is a clear breakpoint where the lowest
%       few singular values are much less than all the others. It is often 
%       good to project these out. Determine the number to be suppressed 
%       and set the variable ns to it.
%    6. Execute recon=aorecon(mvtos,ns) to create the reconstructor.
%    7. Pre multiply recon by the mode-opd & opd-actuator matrices
%    8. Save the workspace. The file can be used as input to WaveTrain 
%       (usually as an argument to a TasatDMModel constructor).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: mod2opd.m 3468 2012-05-17 18:41:01Z keith $

%% BEGIN_CODE

error(nargchk(4,6,nargin));

if nargin==5 && ~isempty(rad)
    annular=true;
    erad=arad/rad;
else
    annular=false;
end
if nargin==6 && ~isempty(ordering)
    % Using individial modes, not orders
    input_modes = true;
    fprintf(1,'mod2opd: using individual modes and specified Zernike ordering\n');
else
    input_modes = false;
    fprintf(1,'mod2opd: using Zernike radial orders\n');
end
%
nloc=length(xopd);
if (length(yopd) ~= nloc); error('mod2opd: length(xopd) ~= length(yopd)'); end;
%
if (sum(orders_modes)==0)
% No modes given - make G an identity matrix so that a tilt included zonal
% reconstructor will be created.
    fprintf(1,'No modes specified. Returning identity matrix.\n');
    fprintf(1,'Generated modal reconstructor maxtrix will be tilt included zonal reconstructor.\n');
    g=eye(nloc);
    return;
end;
%
ms=[];
ls=[];
for i=1:size(orders_modes,2)
    if orders_modes(i) > 0
        if input_modes
            % Get mode from index and ordering
            [l,m]=TwoParameter(orders_modes(i),ordering);
            ms=[ms m];
            ls=[ls l];
        else
            % Get individual zernike modes from radial orders
            m=orders_modes(i);
            for l=-m:2:m
                ms=[ms m];
                ls=[ls l];
            end
        end
    end
end
nmod=size(ms,2);
%
g=zeros(nloc,nmod);
for i=1:nloc
    x=xopd(i); 
    y=yopd(i);
    theta=atan2(y,x);
    r=sqrt(x^2+y^2)/rad;
    for j=1:nmod
        if(~annular)
            %g(i,j)=RealZernike(ls(j),ms(j),r,theta);
            g(i,j)=RealZernike(ls(j),ms(j),r,theta,3); % Use Noll normalization (orthonormal)
        else
            g(i,j)=AnnularZernike(ls(j),ms(j),r,theta,erad);
        end
    end
end
%
% Mirror displacment = 1/2 OPD
g=g/2.0;
%
fprintf(1,'The mode-opd matrix was computed.\n');
