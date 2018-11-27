% out=fftshift(in) : performs the same as fftshift for normal matlab arrays, but here also for dip_image type objects
function out=ifftshift(in)

%safety unequal sizes
if isa(in,'dip_image') 
    out=expanddim(dip_image(ifftshift(double(in))),ndims(in));
else
    out=builtin('ifftshift',in);    % call the builtin version
end
