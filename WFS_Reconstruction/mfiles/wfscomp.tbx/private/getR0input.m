function [s]=getR0input(frameinfo, stats, lambda);
% SYNTAX:
% [s]=getR0input(frameinfo, stats, lambda)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: getR0input.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE
 
if (~missoa(frameinfo))
   ndatasets = length(frameinfo);
   s.xslp = zeros(length(frameinfo(1).slopex),ndatasets);
   s.yslp = s.xslp; 
   for kk=1:1:ndatasets;
      s.xslp(:,kk) = [frameinfo(kk).slopex - stats.slopexAvg] .* (2*pi/lambda);
      s.yslp(:,kk) = [frameinfo(kk).slopey - stats.slopeyAvg] .* (2*pi/lambda);
   end
else
   ndatasets = size(frameinfo.slopex,2);
   s.xslp = frameinfo.slopex;
   s.yslp = frameinfo.slopey;
   for kk=1:1:ndatasets;
      s.xslp(:,kk) = [s.xslp(:,kk) - stats.slopexAvg] .* (2*pi/lambda);
      s.yslp(:,kk) = [s.yslp(:,kk) - stats.slopeyAvg] .* (2*pi/lambda);
   end
end;
