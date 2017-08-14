% [nimg,mykernel]=DampEdgeOutside(img, border, mykernel) : extrapolates the data by filling in blurred information outside the edges. This is a bit like DampEdge but using a normalized convolution with a kernel
% img : image to extrapolate from
% border : percentage of border pixels to add
% mykernel : A kernel can be provided. Default: 1/r^3 kernel of support norm(border size) * sqrt(2)
% usepixels : number of border pixels in the image to use on each edge. Default: 2. Zero means use all pixels
%
% Author: Rainer Heintzmann, IPHT Jena
%
% Example:
% a=readim;
% [damped,kernel]=DampEdgeOutside(a)
% ft(a)
% ft(damped)

function [nimg,mykernel]=DampEdgeOutside(img, border, mykernel, usepixels)
if nargin< 2 || isempty(border)
    border=0.1;
end
if nargin < 4
    usepixels=2;
end
kernelpower=3;
rborder=ceil(border*size(img));
newsize=size(img)+rborder;
if nargin < 3 || isempty(mykernel)
    %mykernel=(1-rr(newsize)/(norm(rborder)/1.0));
    mykernel=(1.0/rr(newsize))-1.0/norm(rborder*sqrt(2));
    mykernel(mykernel < 0) =0;
    mykernel=mykernel^kernelpower;
    mykernel(MidPosX(mykernel),MidPosY(mykernel))=0; % This does not matter as it only applies to points directly having data.
end
transfer=ft(mykernel);
wimg=newim(size(img))+1.0;
nimg=img;
if usepixels > 0
    wimg(usepixels:end-usepixels,usepixels:end-usepixels)=0;
    nimg(usepixels:end-usepixels,usepixels:end-usepixels)=0;
end
nimg=extract(nimg,newsize);
wimg=extract(wimg,newsize); % just mark every pixel

nimg2=real(ift2d(ft2d(nimg) .* transfer));
wimg2=real(ift2d(ft2d(wimg) .* transfer));

nimg(~wimg) = nimg2(~wimg) ./ wimg2(~wimg);
if usepixels > 0
    nimg(floor(rborder/2):end-ceil(rborder/2),floor(rborder/2):end-ceil(rborder/2))=img;
end
