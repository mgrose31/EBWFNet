function X = pinvLM(A,varargin)
% SYNTAX:
% X = pinvLM(A,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PINV   Pseudoinverse.
%   X = PINV(A) produces a matrix X of the same dimensions
%   as A' so that A*X*A = A, X*A*X = X and A*X and X*A
%   are Hermitian. The computation is based on SVD(A) and any
%   singular values less than a tolerance are treated as zero.
%   The default tolerance is MAX(SIZE(A)) * NORM(A) * EPS.
%
%   PINV(A,TOL) uses the tolerance TOL instead of the default.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEE ALSO:  RANK.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
% (c) 1984-2002 The MathWorks, Inc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: pinvLM.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

[m,n] = size(A);
if n > m
   X = pinvLM(A',varargin{:})';
else
   [U,S,V] = svd(A,0);
   if m > 1, s = diag(S);
      elseif m == 1, s = S(1);
      else s = 0;
   end
   clear S;pack
   if nargin == 2
      tol = varargin{1};
   else
      tol = max(m,n) * max(s) * eps;
   end
   clear m n; pack
   r = sum(s > tol);
   clear tol;pack
   if (r == 0)
      X = zeros(size(A'));
   else
      s = diag(ones(r,1)./s(1:r));
      X = V(:,1:r)*s*U(:,1:r)';
   end
end
