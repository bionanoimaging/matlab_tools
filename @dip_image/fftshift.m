% out=fftshift(in) : performs the same as fftshift for normal matlab arrays, but here also for dip_image type objects
function out=fftshift(in)

%safety unequal sizes
if isa(in,'dip_image') 
    out=expanddim(dip_image(fftshift(double(in))),ndims(in));
else
    out=builtin('fftshift',in);    % call the builtin version
end
