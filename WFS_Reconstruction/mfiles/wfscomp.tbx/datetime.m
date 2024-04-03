function datevector = datetime(filename);
% SYNTAX:
% datevector = datetime(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% Input filename must have date and time coded in the final 12 characters,
% as in this example: /ablact/hrwfs/051398/81321133.335 
%   8132 = (199)8, day 132; 
%   1133.335 =  11 hours, 33 minutes, 33 seconds, + 0.5 seconds
% This function returns the matlab date vector using datevec. 
% NOTE: Since the leading digit is presumed to be the last digit of the 
%       year, THIS FUNCTION WILL ONLY WORK correctly THROUGH THE YEAR 2006.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% filename [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS: 
% datavector [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: datetime.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE
 
   filename
   ln = length(filename);
   d1 = ln - 11;
   d2 = ln - 10;
   year = str2num(filename(d1));
   if (year < 6.5)
     year = 2000 + year;
   else
     year = 1990 + year;
   end
   yyyy = num2str(round(year));
   yearnum = datenum(['01-Jan-',yyyy]);
   d2 = ln - 10;
   d3 = d2 + 2;
   day = str2num(filename(d2:d3));
   daynumber = yearnum + day;
   [Y,M,D,h,n,s] = datevec(daynumber);
   t1 = ln - 7;
   t2 = t1 + 1;
   t3 = t1 + 2;
   t4 = t1 + 3;
   t5 = t1 + 5;
   t6 = t1 + 6;
   hh = filename(t1:t2);
   nn = filename(t3:t4);
   ss = filename(t5:t6);
   timeform = [hh,':',nn,':',ss];
   [y,m,d,H,N,S] = datevec(timeform);
   sec = S;
   tenthsec = filename(t1+7);
   dsec = str2num(tenthsec)/10;
   sec = round(sec + dsec);
   S = sec;
   datenumber = datenum(Y,M,D,H,N,S);
   datevector = datevec(datenumber);

