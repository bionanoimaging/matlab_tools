% img=normalizedConvolution(img,mymask,extraBorder,kernel,kernelpower) : a convolution considering a mask. Can be used for inpainting. 
% img : image to process
% mymask : mask in which to keep the image information. Its complement is used for the inpainting
% extraBorder : this can be used to avoid wrap-around effects. If negative, the image will be cut down to original size afterwards.
% kernel : optional argument. Default is a 1/r^4 kernel
% kernelpower : optional argument to modify the kernelstrength, if kernel is left empty []
%
% Example:
% normalizedConvolution(readim,readim > 100,-40)
%
% see also: DampEdgeOutside
%
function img=normalizedConvolution(img,mymask,extraBorder,kernel,kernelpower)
if nargin < 3 || isempty(extraBorder)
    extraBorder=0;
end
if nargin < 5
    kernelpower=4;
end

img=img+0.0;
img(~mymask)=0.0;
origsize=size(img);
newsize=origsize+abs(extraBorder);
if extraBorder ~= 0
    mymask = extract(mymask,newsize);
    img=extract(img,newsize);
end
if nargin < 4 || isempty(mykernel)
    %mykernel=(1-rr(newsize)/(norm(rborder)/1.0));
    mykernel=(1.0/rr(newsize))-1.0/norm(newsize*sqrt(2));
    mykernel(mykernel < 0) =0;
    mykernel=mykernel^kernelpower;
    mykernel(MidPosX(mykernel),MidPosY(mykernel))=1.0; % This does not matter as it only applies to points directly having data.
end
transfer=ft(mykernel);

nimg2=real(ift2d(ft2d(img) .* transfer));
wimg2=real(ift2d(ft2d(mymask) .* transfer));

img(~mymask) = nimg2(~mymask) ./ wimg2(~mymask);

if extraBorder < 0
    img=extract(img,origsize);
end
