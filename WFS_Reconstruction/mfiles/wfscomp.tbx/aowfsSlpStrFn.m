function [dxxx,dxxy,dyyx,dyyy,dxxxy,dyyxy]=aowfsSlpStrFn2(slopex,slopey, ...
    xsub,ysub,xcenter,ycenter,spacing,inradius);
% SYNTAX:
% [dxxx,dxxy,dyyx,dyyy,dxxxy,dyyxy]=aowfsSlpStrFn2(slopex,slopey,xsub, ...
%                                   ysub, xcenter,ycenter,spacing,inradius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% IN order to locate the correct center radius need to perform the 
% following conversion; compute dxxx and dyyx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% slopex [ ] =
% slopey [ ] = 
% xsub [ ] = 
% ysub [ ] = 
% xcenter [ ] = 
% ycenter [ ] = 
% spacing [ ] = 
% inradius [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% dxxx [ ] =
% dxxy [ ] = 
% dyyx [ ] = 
% dyyy [ ] = 
% dxxxy [ ] = 
% dyyxy [ ] = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: aowfsSlpStrFn.m 3063 2010-10-08 20:42:07Z amoran $

%% BEGIN_CODE
 
aa = find(xsub > (xcenter-spacing) & xsub < (xcenter+spacing));
tmp = ysub(aa);
aa1 = find((ysub(aa)-ycenter) > 0);
aa2 = find((ysub(aa)-ycenter) < 0);
aax1 = min(tmp(aa1));
aax2 = max(tmp(aa2));

if (isempty(aax1) | isempty(aax2))
  dxxx = [];
  dxxy = [];
  dyyx = [];
  dyyy = [];
  dxxxy = [];
  dyyxy = [];
  return
end

tmpx1 = find((ysub < (aax1+spacing/2)) & (ysub > (aax1-spacing/2)));
tmpx2 = find((ysub < (aax2+spacing/2)) & (ysub > (aax2-spacing/2)));

nfms = size(slopex,2);
nsub = min(length(tmpx1),length(tmpx2));
dxxx = zeros(nsub,nfms);
dyyx = zeros(nsub,nfms);
for kkk=1:2
   if kkk==1
      tmpx=tmpx1;
   else
      tmpx=tmpx2;
   end
   slopext=slopex(tmpx,:);
   slopeyt=slopey(tmpx,:);
   for kk=1:nsub
       dxxxt = sum((slopext(kk:nsub,:)-slopext(1:(nsub-kk+1),:)).^2,1);
       dyyxt = sum((slopeyt(kk:nsub,:)-slopeyt(1:(nsub-kk+1),:)).^2,1);
       dxxx(kk,:) = dxxx(kk,:) + dxxxt(1,:)./(nsub-kk+1);
       dyyx(kk,:) = dyyx(kk,:) + dyyxt(1,:)./(nsub-kk+1);
   end
end
dxxx = dxxx/2;
dyyx = dyyx/2;

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

nfms=size(slopex,2);
nsub=min(length(tmpy1),length(tmpy2));
dxxy=zeros(nsub,nfms);
dyyy=zeros(nsub,nfms);

for kkk=1:2
   if kkk==1
      tmpy=tmpy1;
   else
      tmpy=tmpy2;
   end;
   slopext=slopex(tmpy,:);
   slopeyt=slopey(tmpy,:);
   for kk=1:nsub
      dxxyt = sum((slopext(kk:nsub,:)-slopext(1:(nsub-kk+1),:)).^2,1);
      dyyyt = sum((slopeyt(kk:nsub,:)-slopeyt(1:(nsub-kk+1),:)).^2,1);
      dxxy(kk,:) = dxxy(kk,:) + dxxyt(1,:)./(nsub-kk+1);
      dyyy(kk,:) = dyyy(kk,:) + dyyyt(1,:)./(nsub-kk+1);
   end;
end;
dxxy = dxxy/2;
dyyy = dyyy/2;

% 
% compute dxxxy and dyyxy
%

b1=-1*(aay2-aax1);
b2=-1*(aay1-aax2);

tmpxy1 = find(xsub+b1>=ysub-spacing/2 & xsub+b1<=ysub+spacing/2);
tmpxy2 = find(xsub+b2>=ysub-spacing/2 & xsub+b2<=ysub+spacing/2);

nfms = size(slopex,2);
nsub = min(length(tmpxy1),length(tmpxy2));
dxxxy = zeros(nsub,nfms);
dyyxy = zeros(nsub,nfms);

for kkk=1:2
   if kkk==1
      tmpxy=tmpxy1;
   else
      tmpxy=tmpxy2;
   end
   slopext=slopex(tmpxy,:);
   slopeyt=slopey(tmpxy,:);
   for kk=1:nsub
      dxxxyt = sum((slopext(kk:nsub,:)-slopext(1:(nsub-kk+1),:)).^2,1);
      dyyxyt = sum((slopeyt(kk:nsub,:)-slopeyt(1:(nsub-kk+1),:)).^2,1);
      dxxxy(kk,:) = dxxxy(kk,:) + dxxxyt(1,:)./(nsub-kk+1);
      dyyxy(kk,:) = dyyxy(kk,:) + dyyxyt(1,:)./(nsub-kk+1);
   end
end
dxxxy = dxxxy/2;
dyyxy = dyyxy/2;
