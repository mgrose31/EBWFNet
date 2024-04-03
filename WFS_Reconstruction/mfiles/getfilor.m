function [fileout] = getfilor(dirname, cond)
% SYNTAX:
% [fileout] = getfilnm(dirname, cond)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% determines all files in directory dirname such that a portion of their 
% name is at least one of the strings contained in the cell cond. -- this 
% is in contrast to getfilnm which determines all filenames which match 
% all of the criteria set in cond.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% dirname [ ] =
% cond [ ] =
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% fileout [ ] =
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: getfilor.m 3051 2010-10-01 20:33:26Z amoran $

%% BEGIN_CODE

warning off;

direct = dir(dirname);
[n1,n2] = size(direct);
[n3,n4] = size(cond);

for i = 1:n1
  for j = 1:n4
    lstr = length(cond{j});
    lname = length(direct(i).name);
    if (strcmp(cond{j}(1),'*') == 1)
      lpos = findstr(cond{j}(2:lstr), direct(i).name);
      if (lpos == lname-lstr+2) 
        strpos{i,j} = lpos;
      else
        strpos{i,j} = [];
      end
    elseif (strcmp(cond{j}(lstr),'*') == 1)
      lpos = findstr(cond{j}(1:lstr-1), direct(i).name);
      if (lpos == 1) 
        strpos{i,j} = lpos;
      else
        strpos{i,j} = [];
      end
    else
      strpos{i,j} = findstr(cond{j}, direct(i).name);
    end 
  end
end

findit = ones(n1,1);
for i = 1:n1
  tmp = 0;
  for j = 1:n4
%
%    findit(i) = ~isempty(strpos{i,j})*findit(i);
%
     tmp = max(~isempty(strpos{i,j}),tmp);
  end
  findit(i)=tmp;
end

count = 1
clear fileout;
for i = 1:n1
  if (findit(i) == 1)
    fileout{count} = direct(i).name;
    count = count + 1;
  end
end

warning on;
