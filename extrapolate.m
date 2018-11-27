% res=extrapolate(img,amask): interpolates (extrapolates) values for the locations not in the amask image. Is also a bit like hole filling. Sometimes also called inpainting. 
% img: input image
% amask : mask image
%
% Example:
% extrapolate(readim,readim<100)

function res=extrapolate(img,amask)

myX=xx(img,'corner');myY=yy(img,'corner');

Vdat=double(img(amask));
Xdat=double(myX(amask));
Ydat=double(myY(amask));

[Xq,Yq] = meshgrid(0:size(img,1)-1,0:size(img,2)-1);

res=mat2im(griddata(Xdat,Ydat,Vdat,Xq,Yq,'cubic'));
