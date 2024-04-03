function dmgtname(dostring)
% SYNTAX:
% dmgtname(dostring)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% dostring [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: dmgtname.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE

persistent prevdir;
curdir = cd;
if (~isempty(prevdir))
   if (~exist(prevdir,'dir'))
      prevdir = [];
   end;
end;
if (isempty(prevdir))
   wtdir = [getenv('WT_DIR'),'\..\predata'];
   if (exist(wtdir, 'dir'))
      cd(wtdir);
   else
      wtdir = [getenv('WT_DIR'),'\predata'];
      if (exist(wtdir, 'dir'))
         cd(wtdir);
      end;
   end;
else
   cd(prevdir);
end;
if (strcmp(dostring,'load'))
   [Filename,Pathname] = uigetfile('*.mat','Load AO Configuration...');
else
   [Filename,Pathname] = uiputfile('*.mat','Save AO Configuration...');
end;
if (~Pathname)
   Pathname = [cd,'\'];   
end;
if (Filename)
   prevdir = Pathname;
   if (strcmp(dostring,'load'))
      dmload([Pathname, Filename]);
      dmplot(1);
   else
      dmsave([Pathname, Filename]);
   end;
end;
cd(curdir);
