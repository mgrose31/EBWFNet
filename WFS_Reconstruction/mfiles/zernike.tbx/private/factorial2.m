function [f]=factorial2(in)
% SYNTAX:
% [f]=factorial2(in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: factorial2.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE 

if (in<=0) 
    f=1;
    return; 
end;
if (in>10)
    f=factorial(in);
    return;
else
    switch in
        case 1
            f=1; return;
        case 2
            f=2; return;
        case 3
            f=6; return;
        case 4
            f=24; return;
        case 5
            f=120; return;
        case 6
            f=720; return;
        case 7
            f=5040; return;
        case 8
            f=40320; return;
        case 9
            f=362880; return;
        case 10
            f=3628800; return;
    end;
end
return;