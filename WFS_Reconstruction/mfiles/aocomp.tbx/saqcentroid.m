function [centroid] = saqcentroid(sa)
% SYNTAX: 
% [centroid] = saqentroid(sa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS: 
% sa [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS: 
% centroid [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: saqcentroid.m 3061 2010-10-07 21:13:39Z amoran $

%% BEGIN_CODE

% avg1 = sum(sa)  ./ size(sa,1);
% avg2 = sum(sa.')./ size(sa,2);
% num1 = sum(avg1 .* [1:size(sa,2)]);
% num2 = sum(avg2 .* [1:size(sa,1)]);
% den1 = sum(avg1);
% den2 = sum(avg2);
% if ((den1*den2) == 0)
%    centroid = [mean([1,size(sa,1)]), mean([1,size(sa,2)])];
% else
%    centroid = [num2./den2, num1./den1];
% end;


sizesa = size(sa);
avg1 = squeeze(shiftdim(mean(sa,1),-1));
avg2 = squeeze(mean(sa,2));
x = 1:sizesa(2);
y = 1:sizesa(1);
num1 = x * avg1;
num2 = y * avg2;
den1 = sum(avg1,1);
den2 = sum(avg2,1);
idx = den1.*den2 == 0;
centroid = 0*[num1; num2];
centroid(1,idx) = (1+sizesa(1))/2;
centroid(2,idx) = (1+sizesa(2))/2;
idx = ~idx;
centroid(:,idx) = [num2(idx)./den2(idx); num1(idx)./den1(idx)];
centroid = centroid.';