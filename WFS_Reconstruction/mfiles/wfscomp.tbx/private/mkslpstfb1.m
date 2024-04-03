function [dxxx,dxxy,dyyx,dyyy,dxxxy,dyyxy]=mkslpstfb1(slopex,slopey,xsub,ysub,xcenter,ycenter,spacing,inradius);
% SYNTAX:
% [dxxx,dxxy,dyyx,dyyy,dxxxy,dyyxy]=mkslpstfb1(slopex,slopey,xsub,ysub, ...
%                                   xcen ter,ycenter,spacing,inradius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% IN order to locate the correct center radius need to perform the 
% following conversion -- tested using: compute dxxx and dyyx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: mkslpstfb1.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE
 
aa = find(xsub>xcenter-spacing & xsub<xcenter+spacing);
tmp=ysub(aa);
aa1 = find(ysub(aa)-ycenter>0);
aa2 = find(ysub(aa)-ycenter<0);
aax1 = min(tmp(aa1));
aax2 = max(tmp(aa2));

if(isempty(aax1)|isempty(aax2))
  dxxx=0;
  dxxy=0;
  dyyx=0;
  dyyy=0;
  dxxxy=0;
  dyyxy=0;
  return
end

tmpx1 = find(ysub<aax1+spacing/2 & ysub>aax1-spacing/2);
tmpx2 = find(ysub<aax2+spacing/2 & ysub>aax2-spacing/2);

nsub=min(length(tmpx1),length(tmpx2));
nfms=size(slopex,2);
dxxx=zeros(1,nsub);
dyyx=zeros(1,nsub);
for kk=1:1:2
if kk==1
   tmpx=tmpx1;
else
   tmpx=tmpx2;
end

for jj=1:1:nfms
   slopext=slopex(tmpx,jj);
   slopeyt=slopey(tmpx,jj);
   for kk=1:1:nsub
       kk=kk-1;
       dxxxt=0;
       dyyxt=0;
       cnt=0;
       for nn=1:1:nsub-kk
           dxxxt=(slopext(nn+kk)-slopext(nn))^2+dxxxt;
           dyyxt=(slopeyt(nn+kk)-slopeyt(nn))^2+dyyxt;
           cnt=cnt+1;
       end
       dxxx(kk+1)=dxxxt/cnt+dxxx(kk+1);
       dyyx(kk+1)=dyyxt/cnt+dyyx(kk+1);
   end
end
end
dxxx=dxxx/(nfms*2);
dyyx=dyyx/(nfms*2);


%
% compute dxxy and dyyy
%
aa = find(ysub>ycenter-spacing & ysub<ycenter+spacing);
tmp=xsub(aa);
aa1 = find(xsub(aa)-xcenter>0);
aa2 = find(xsub(aa)-xcenter<0);
aay1 = min(tmp(aa1));
aay2 = max(tmp(aa2));

tmpy1 = find(xsub<aay1+spacing/2 & xsub>aay1-spacing/2);
tmpy2 = find(xsub<aay2+spacing/2 & xsub>aay2-spacing/2);

nsub=min(length(tmpy1),length(tmpy2));
nfms=size(slopex,2);
dxxy=zeros(1,nsub);
dyyy=zeros(1,nsub);

for kk=1:1:2
if kk==1
   tmpy=tmpy1;
else
   tmpy=tmpy2;
end

for jj=1:1:nfms
   slopext=slopex(tmpy,jj);
   slopeyt=slopey(tmpy,jj);
   for kk=1:1:nsub
       kk=kk-1;
       dxxyt=0;
       dyyyt=0;
       cnt=0;
       for nn=1:1:nsub-kk
           dxxyt=(slopext(nn+kk)-slopext(nn))^2+dxxyt;
           dyyyt=(slopeyt(nn+kk)-slopeyt(nn))^2+dyyyt;
           cnt=cnt+1;
       end
       dxxy(kk+1)=dxxyt/cnt+dxxy(kk+1);
       dyyy(kk+1)=dyyyt/cnt+dyyy(kk+1); 
   end
end
end
dxxy=dxxy/(nfms*2);
dyyy=dyyy/(nfms*2);

% 
% compute dxxxy and dyyxy
%

b1=-1*(aay2-aax1);
b2=-1*(aay1-aax2);

tmpxy1 = find(xsub+b1>=ysub-spacing/2 & xsub+b1<=ysub+spacing/2);
tmpxy2 = find(xsub+b2>=ysub-spacing/2 & xsub+b2<=ysub+spacing/2);

nsub=min(length(tmpxy1),length(tmpxy2));
nfms=size(slopex,2);
dxxxy=zeros(1,nsub);
dyyxy=zeros(1,nsub);

for kk=1:1:2
if kk==1
   tmpxy=tmpxy1;
else
   tmpxy=tmpxy2;
end

for jj=1:1:nfms
   slopext=slopex(tmpxy,jj);
   slopeyt=slopey(tmpxy,jj);
   for kk=1:1:nsub
       kk=kk-1;
       dxxxyt=0;
       dyyxyt=0;
       cnt=0;
       for nn=1:1:nsub-kk
           dxxxyt=(slopext(nn+kk)-slopext(nn))^2+dxxxyt;
           dyyxyt=(slopeyt(nn+kk)-slopeyt(nn))^2+dyyxyt;
           cnt=cnt+1;
       end
       dxxxy(kk+1)=dxxxyt/cnt+dxxxy(kk+1);
       dyyxy(kk+1)=dyyxyt/cnt+dyyxy(kk+1);
   end
end
end
dxxxy=dxxxy/(nfms*2);
dyyxy=dyyxy/(nfms*2);


