function vmat=slope2mod(orders_modes,xsub,ysub,rad,tol,ordering)
% SYNTAX: 
% vmat=slope2mod(orders,xsub,ysub,rad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: 
% slope2mod computes the slope to Zernike 
% Coefficient matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% orders [ ] = radial orders of Zernike modes to be fit. All angular 
%              orders are fit for a given radial order.
% modes [] = The nmode x 1 vector of Zernike modes to be used. Requires the
%              'ordering' parameter to be set
% xsub,ysub [ ] = Locations of the subapetures in output space.
% rad [ ] = Radius of the DM / beam/ aperture
% tol [ ] =  Least squares / orthononality tolerance.
%       >=1 ==> overlap integral
%       <=0 ==> least squares
%       1>tol>0 ==> test matrix
% ordering = Zernike ordering scheme (optional); 1=Noll, 2=Malacara, 3=Wyant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% vmat [matrix] = The block matrix formed as follows, but then re-arranged 

% /                                                                          \
% |   v^1_x(R_1)  ...   v^1_x(r_nsub)  |   v^1_y(R_1)  ...    v^1_y(R_nsub)  |
% |       :       :::        :         |       :       :::        :          |
% | v^nmod_x(R_1) ... v^nmod_x(r_nsub) | v^nmod_y(R_1) ...  v^nmod_y(R_nsub) |
% \                                                                          /
% where:
% v^i_x(R_j) is the value of the X component of the Gavrielides
% vector polynomial corresponding to the ith Zernike mode at the jth
% subaperture location
% v^i_x(R_j) = (dx/dr, dx/dtheta)_(R_j) . (Gr(i), Gt(i))*
%
% v^i_y(R_j) is the value of the Y component of the Gavrielides
% vector polynomial corresponding to the ith Zernike mode at the jth
% subaperture location
%   v^i_y(R_j) = (dy/dr, dy/dtheta)_(R_j) . (Gr(i), Gt(i))*
%
% Notes:
% 1i is not set - it takes its default value of sqrt(-1);
% We first calculate the above matrix that will give coefficients of 
% general (complex) Zernikes. We then modify this matrix to give 
% coefficients of real Zernikes such that
%
%           /    a^C_m,l+a^C_m,-l  l<0  cos(theta) real term
% a^R_m,l = |    a^C_m,l           l=0 (spherically symmetric) should be real anyway
%           \ i*(a^C_m,l-a^C_m,-l) l>0  sin(theta) real term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: slope2mod.m 3468 2012-05-17 18:41:01Z keith $

%% BEGIN_CODE

error(nargchk(4,6,nargin))

if nargin < 5 || isempty (tol)
    tol = 0.1;
end
if nargin==6 && ~isempty(ordering)
    % Using individial modes, not orders
    input_modes = true;
    fprintf(1,'slope2mod: using individual modes and specified Zernike ordering\n');
else
    input_modes = false;
    fprintf(1,'slope2mod: using Zernike radial orders\n');
end
%
nsub=size(xsub,1);
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
nmode_in=size(ms,2); % Number of input modes
% Need symmetric modes so we can convert from complex to real
nls = zeros(size(ls));
for i=1:nmode_in
    l = ls(i);
    m = ms(i);
    for j=i:nmode_in
        if ms(j) == m && ls(j) == -l
            nls(i) = j;
            nls(j) = i;
            break
        end
    end
    if nls(i) == 0
        ms = [ms m];
        ls = [ls -l];
        nls(i) = length(ls);
    end
end
nmode=size(ms,2); % Number of modes after additions for symmetry   
%
% Create vmat, complex Gavrieldes matrix, i=1:nmodes, j=1:nsubaps
% [(dx/dr, dx/dtheta)_(R_j) . (Gr(i), Gt(i))* | (dy/dr, dy/dtheta)_(R_j) . (Gr(i), Gt(i))*]
%
vmat = nan*ones(nmode,2*nsub);
for imode=1:nmode
    m=ms(imode);
    l=ls(imode);
    %c = sqrt(2*(m+1)); % Use Noll normalization (orthonormal) for complex Zernikes
    c = sqrt((m+1)/2);
    for jsub=1:nsub
        x=xsub(jsub)/rad;
        y=ysub(jsub)/rad;
        r=sqrt(x^2+y^2);
        theta=atan2(y,x);
        % Uses Liyang's Matlab scripts
        Gr = c*KradGZ(abs(l),m,r)*exp(1i*l*theta);
        Gt = c*1i*l*SigmaGZ(abs(l),m,r)*exp(1i*l*theta)/r;
        % Conjugate & convert to Cartesian
        % conj(Gr)dx/dr + conj(Gt)dx/dt
        Gx = rad*(conj(Gr)*cos(theta)-conj(Gt)*y);
        % conj(Gr)dy/dr + conj(Gt)dy/dt
        Gy = rad*(conj(Gr)*sin(theta)+conj(Gt)*x);
        vmat(imode,jsub)=Gx;
        vmat(imode,jsub+nsub)=Gy;
    end
end
%
% Normalize
vmat = vmat/nsub;
%
% Orthogonality check (coarse mesh, obscuration, etc) and least squares fix.
% Matrix:
% [Sum{ (dZ(k)/dR)_(R_j) . (Gr(i), Gt(i))* }j=1:nsubaps], i,k=1:nmodes
clear m l r t x y theta
zmat = nan*ones(nmode,nmode);
for imode=1:nmode % Loop over Zernike modes for Gavrielides polynomials
    m_i = ms(imode);
    l_i = ls(imode);
    c_i = sqrt(m_i+1);
    
    for kmode = 1:nmode % Loop over Zernike modes for Zernike derivatives
        m_k = ms(kmode);
        l_k = ls(kmode);
        c_k = sqrt(m_k+1); % Use Noll normalization (orthonormal) for complex Zernikes
        if l_k == 0
            c_k=c_k/sqrt(2);
        end
        %
        zmat(imode,kmode) = complex(0.0);
        for jsub=1:nsub % Loop over sub apertures
            x=xsub(jsub)/rad;
            y=ysub(jsub)/rad;
            r=sqrt(x^2+y^2);
            theta=atan2(y,x);
            % Uses Liyang's Matlab scripts
            Gr = c_i*KradGZ(abs(l_i),m_i,r)*exp(1i*l_i*theta);
            Gt = c_i*1i*l_i*SigmaGZ(abs(l_i),m_i,r)*exp(1i*l_i*theta)/r;
            % Keith added Zernike derivatives
            Zr = c_k*dZdr(abs(l_k),m_k,r)*exp(1i*l_k*theta);
            Zt = c_k*1i*l_k*RGZ(abs(l_k),m_k,r)*exp(1i*l_k*theta);
            %
            zmat(imode,kmode) = zmat(imode,kmode) + ( Zr*conj(Gr) + Zt*conj(Gt));
        end
    end
end
% Normalize
zmat = zmat/nsub;
% Matrix should be the identity if Zernike derivatives are orthongonal to
% Gavrielides polyniomials & matrix should always be real
if max(abs(imag(zmat(:))))/mean(abs(real(zmat(:)))) > 1e-6
    error('slope2mod.m: Least-squares matrix is complex!')
else
    zmat = real(zmat); % Force to be real
    figure;imagesc(zmat);axis image;colorbar
    title('Gavrielides Least-Squares Matrix (should be identity if derivatives are orthogonal)')
    figure;imagesc(inv(zmat));axis image;colorbar
    title('Inverse Gavrielides Least-Squares Matrix')
    if tol <= 0
        leastSquares = true;
    elseif tol >= 1
        leastSquares = false;
    else
        zerr = zmat-eye(nmode); % Difference from identity
        leastSquares = (max(abs(zerr(:))) > tol);
        clear zerr
    end
end
if leastSquares
    disp('Calculating least-squares solution')
    vmat = zmat\vmat; % inv(zmat)* vmat
end
%
% The above matrix will give coefficients of general (complex) Zernikes
% Convert to give coefficients of real Zernikes
% Create vmatr, real Gavrieldes matrix, i=1:nmodes, j=1:nsubaps
% [(dx/dr, dx/dtheta)_(R_j) . (Gr(i), Gt(i))* | (dy/dr, dy/dtheta)_(R_j) . (Gr(i), Gt(i))*]
vmatr = nan*ones(nmode_in,2*nsub);
for imode=1:nmode_in
    l=ls(imode);
    if(l==0)
        vmatr(imode,1:2*nsub)=vmat(imode,:);
    else
        imodeml = nls(imode);
        if(l<0)
            vmatr(imode,1:2*nsub)=   (vmat(imode,:)+vmat(imodeml,:));
        elseif(l>0)
            vmatr(imode,1:2*nsub)=1i*(vmat(imode,:)-vmat(imodeml,:));
        end
    end
end
%
if max(abs(imag(vmatr(:))))/mean(abs(real(vmatr(:)))) > 1e-6
    disp('Warning: some components of slope to mode matrix are complex')
    disp(['Maximum magnitude of imaginary components is ',num2str(max(abs(imag(vmatr(:)))))])
    disp(['Maximum magnitude of real components is ',num2str(max(abs(real(vmatr(:)))))])
    disp('Matrix forced to be real !!')
end
% vmat=real(vmatr)/2; %Factor of two for OPD -> Mirror displacment
vmat=real(vmatr);
%
fprintf(1,'The slope-mode matrix was computed.\n');
