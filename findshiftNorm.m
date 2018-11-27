% [shiftvec,shifted]=findshiftNorm(img,ref,minmap,doPhaseOnly,doCOM) : Finds the shift between images, but corrects for the relative overlap area
% img : image to align
% ref : reference image to align to
% minmap :  (default 0.5)
% doPhaseOnly : flag, if true, the phase only correlation is calculated
% doCOM : Use the center of mass of the surrounding 5x5 region to determine the subpixel shift
%
% Example:
% img=readim('orka');
% a=extract(img,[100 100]);
% as=extract(shift(img,[20.2 22.7]),[100 100]);
% [shiftvec,shifted]=findshiftNorm(as-mean(as),a-mean(a),0.5,1,1)
% cat(3,a,shifted)

function [shiftvec,shifted]=findshiftNorm(img,ref,minmap,doPhaseOnly,doCOM)
if nargin < 5
    doCOM=0;
end
if nargin < 4
    doPhaseOnly=0;
end
if nargin < 3
    minmap=0.5;
end
mypow=0.5;

fimg=ft(img);
fref=ft(ref);

if doPhaseOnly
    mycor=abs(ift(exp(i*(phase(fimg) -phase(fref)))));
else
    mycor=abs(ift(fimg .* conj(fref)));
end
%ref2=extract(ref.*ref,floor(size(ref)/2)-1,floor(size(ref)/2));  % energy (not summed, when perfectly aligned)
%ref3=extract(ref2,floor(size(ref)/2));  % energy (not summed, when perfectly aligned)
%ref2=ref.*ref;  % energy (not summed, when perfectly aligned)
%ref3=extract(ref2,size(ref)+1,floor(size(ref)/2-1));  % energy (not summed, when perfectly aligned)
ref3=ref.*ref;  % energy (not summed, when perfectly aligned)
ref2=extract(ref3,size(ref)+1);  % energy (not summed, when perfectly aligned)
tsum=sum(ref2);
ref2=ref2/tsum;
ref3=ref3/tsum;
% now a map needs to be created with the lost energy when shifting out of the image by dx,dy
sumx=dip_image(cumsum(double(ref3),2));  % shifting in -x
summx=flipdim(dip_image(cumsum(double(flipdim(ref2,1)),2)),1);  % shifting in +x
sumy=dip_image(cumsum(double(ref3),1));  % shifting in -y
summy=flipdim(dip_image(cumsum(double(flipdim(ref2,2)),1)),2);  %shifting in y
sumxy=dip_image(cumsum(double(sumy),2));  % shifting in -y
summxmy=flipdim(dip_image(cumsum(double(flipdim(summy,1)),2)),1);  % shifting in -y
sumxmy=dip_image(cumsum(double(summy),2));  % shifting in -y
summxy=dip_image(cumsum(double(summx),1));  % shifting in -y


sumtotalmm=flipdim(flipdim(sum(sumx,[],2)+sum(sumy,[],1)-sumxy,1), 2);  % this energy is lost with each -X, -Y shift
sumtotalpp=flipdim(flipdim(sum(summx,[],2)+sum(summy,[],1)-summxmy,1), 2);  % this energy is lost with each -X, -Y shift
sumtotalmp=flipdim(flipdim(sum(sumx,[],2)+sum(summy,[],1)-extract(sumxmy,[size(sumx,1) size(summy,2)]),1),2);  % this energy is lost with each -X, -Y shift
sumtotalpm=flipdim(flipdim(sum(summx,[],2)+sum(sumy,[],1)-extract(summxy,[size(summx,1) size(sumy,2)]),1),2);  % this energy is lost with each -X, -Y shift



totalmap=extract(cat(2,cat(1,sumtotalmm,sumtotalpm),cat(1,sumtotalmp,sumtotalpp)),size(ref),size(ref));

totalmap=1-totalmap;
%totalmap=1-totalmap/sum(ref2);
totalmap(totalmap<minmap)=minmap;  % regularise it to avoid dividing by zero
% assuming the energy to be equally 

if (~doPhaseOnly)
    [mv,shiftvec]=max(mycor./(totalmap.^mypow));
else
    [mv,shiftvec]=max(mycor);
end
if doCOM
    ROI=extract(mycor,[5 5],shiftvec);
    ROI=ROI-min(ROI);
    COMX=sum(xx(ROI) .* ROI)/sum(ROI);
    COMY=sum(yy(ROI) .* ROI)/sum(ROI);
    shiftvec=shiftvec+[COMX COMY];
end

shiftvec = transpose(shiftvec - floor(size(fimg)/2));
if nargout > 1
    shifted = shiftClipFlip(img,-shiftvec);
end

