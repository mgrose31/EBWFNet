function dmhelp
%%if (exist([getenv('WT_DIR'),'\doc\aogeom\aogeom.htm'],'file'))
%%   web([getenv('WT_DIR'),'\doc\aogeom\aogeom.htm']);
if (exist([getenv('WT_DIR')  '\..\..\doc\aogeom\index.htm'],'file'))
   web([getenv('WT_DIR')  '\..\..\doc\aogeom\index.htm#overaogeom'], '-browser');
else
   web('http://www.mza.com/doc/wavetrain/aogeom/main.htm', '-browser');
end;
