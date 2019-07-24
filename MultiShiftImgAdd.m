% [Img,mygrad]=MultiShiftImgAdd(Obj,PosList): places an object (Obj) at a number of positions (PosList).
% Obj: object to be shifted to the positions
% PosList: a list of positions (as a Matrix) at which the object should be placed
%
% Example:
% BdsCoords=rand(3,10); % relative to image size
% AnImg=exp(-(xx.^2+yy.^2)/50)
% MultiShiftImgAdd(AnImg,BdsCoords)
function Img=MultiShiftImgAdd(Obj,PosList) 
% create list of shifts
kx=reshape(dip_image(PosList(1,:)),[1 1 size(PosList,2)]);
ky=reshape(dip_image(PosList(2,:)),[1 1 size(PosList,2)]);
brightness=reshape(dip_image(PosList(3,:)),[1 1 size(PosList,2)]);
allshifts=sum(brightness.*exp(1i*2*pi*(xx(size2d(Obj)).*kx+yy(size2d(Obj)).*ky)),[],3);

Img = real(ift(ft(Obj).* allshifts));
