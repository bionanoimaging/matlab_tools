% [res,midpos]=MidVal(anImg,val,whichdims)  : returns the value in the center of the image (according to the Fourier-conventions where the middle pixel is)
% anImg: input image
% val: value(s) to assign
% whichdims: (optional) which dims the midpoint shall be extracted
% res: middle value (as a double)
% midpos: index where the midpoint was determined (DipImage convention starting with zero)
% 
% Example:
% sum(readim)
% MidVal(ft(readim)*sqrt(prod(size(readim))))

function img=MidValAsg(img,val,whichdims)
if nargin < 3
    whichdims=[1:ndims(img)];
end
idx=cell(1,ndims(img));
for n=1:ndims(img)
    if any(n==whichdims)
        idx{n}=floor(size(img,n)/2);
    else
        idx{n}=':';
    end
end
S=struct('type','()','subs',{idx});
img=subsasgn(img,S,val);
