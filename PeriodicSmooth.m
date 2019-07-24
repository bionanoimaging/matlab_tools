% function p = PeriodicSmooth(img): generate the periodic component of a 2D image which dim the border of an image and avoid discontinuities in Fourier space according to Moisan et al.
%
% img : Input image to smoothen
% eg: 
% Erica=PeriodSmooth(readim());
% p=ft(Erica);
% diff=p-ft(readim());
% diffRealSpace=ift(diff);
%
% Ref.: Moisan, J Math Imaging Vis (2011) 39: 161–179, DOI 10.1007/s10851-010-0227-1

function p = PeriodicSmooth(img)
sz=size(img);
mysum=newim(sz);
den=-2*ndims(img);
for d=1:ndims(img)
    low=SubSlice(img,d,0);
    high=SubSlice(img,d,sz(d)-1);
    mysum=SubSliceAsg(mysum,d,0,high-low + SubSlice(mysum,d,0));
    mysum=SubSliceAsg(mysum,d,sz(d)-1,low-high + SubSlice(mysum,d,sz(d)-1));
    den=den+2*cos(2*pi*ramp(dimVec(d,sz(d),length(sz)),d,'freq'));
end
den=ft(mysum)./den;
den=MidValAsg(den,0);
den=ift(den);
p=img-den;
