% [res,midpos]=MidVal(anImg)  : returns the value in the center of the image (according to the Fourier-conventions where the middle pixel is)
% anImg: input image
% res: middle value (as a double)
% midpos: index where the midpoint was determined (DipImage convention starting with zero)
% 
% Example:
% sum(readim)
% MidVal(ft(readim)*sqrt(prod(size(readim))))

function img=MidValAsg(img,val)
idx=cell(1,ndims(img));
for n=1:ndims(img)
    idx{n}=floor(size(img,n)/2);
end
S=struct('type','()','subs',{idx});
img=subsasgn(img,S,val);
