% res=fftshift2d(animg,varargin): applied fftshift to each 2d slice in a 3d array.
% animg: needs to be 3d
%
function res=fftshift2d(animg,varargin)
if ndims(animg) > 2
    if isa(animg,'dip_image') || (isa(animg,'cuda') && isDip(animg))
        off=1;
    else
        off=0;
    end
    res={};
    for n=1:size(animg,3)
        res{n}=fftshift(animg(:,:,n-off));
    end
    res=cat(3,res{:});
else
    res=fftshift(animg);
end
