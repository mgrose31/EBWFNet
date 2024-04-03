function [l,n]=Old_Noll(index)
% SYNTAX:
% [l,n]=Old_Noll(index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% This used to be called 'Noll' but is totally different to Noll's ordering
% Now considered depreciated and should not be used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% index [ ] = 1D Zernike index (starting with tilts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% n [ ] = radial order (starting at 1 for tilts)
% l [ ] = angular order (one of -n:2:n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: Old_Noll.m 3036 2010-09-23 21:22:50Z amoran $

%% BEGIN_CODE 

l=0;n=0;
cnt=0;
for ii=0:index;
    for jj=-ii:2:ii;
        if (cnt==index)
            l=jj; n=ii;
            return;
        end
        cnt=cnt+1;
    end;
end
return;