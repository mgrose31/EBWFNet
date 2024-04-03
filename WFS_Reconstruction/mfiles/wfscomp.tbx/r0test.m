function [r,sig,ahat]=r0test(s,rad,w,c,guess,iter);
% SYNTAX:
% [r,sig,ahat]=r0test(s,rad,w,c,guess,iter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% s [ ] = slopesarray 1-d array containing frames of x-slope
%         followed by y-slope values.
% num [ ] = the number of frames in the array slopearray
% w [ ] = linear transformation mapping zernikes to slopes
% c [ ] = the variance matrix of zernikes
% rad [ ] = the radius in meters of the pupil aperture
% guess & iter [ ] = self explanatory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% r [ ] = 
% sig [ ] = 
% ahat [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: r0test.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE
 
d = w*c*w.';
nn = size(d,1);
num=size(s.xslp,2);
for ii = 1:1:iter
   h = c*w.'*(((d + guess*eye(size(d,1))))^-1);
   vara=0;
   vars=0;
   for kk=1:1:num
      slopes = [s.xslp(:,kk).', s.yslp(:,kk).'];
      vars = vars + sum(slopes.^2);
      numsubap = 2*(length(slopes)/2);
      ahat{kk} = h*(slopes.');
      vara=vara+sum(ahat{kk}.^2);
   end
   vars=vars/num;
   vara=vara/num;
   a = [trace(h*d*h.') trace(h*h.');trace(d) numsubap];
   b = [vara vars].';
   x=pinv(a)*b; 
   guess = x(2)/x(1);
   det(a);
   cond(a);
   r=2.*rad/(x(1)^(3/5));
   sig=x(2);
end

