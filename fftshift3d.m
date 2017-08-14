function res=fftshift3d(animg,varargin)
if ndims(animg) > 3
    if isa(animg,'dip_image')
        off=1;
    else
        off=0;
    end
    res={};
    for n=1:size(animg,4)
        res{n}=fftshift(animg(:,:,:,n-off));
    end
    res=cat(4,res{:});
else
    res=fftshift(animg);
end
