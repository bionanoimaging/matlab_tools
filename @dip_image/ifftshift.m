% out=fftshift(in) : performs the same as fftshift for normal matlab arrays, but here also for dip_image type objects
function out=ifftshift(in)

%safety unequal sizes
if isa(in,'dip_image') 
    out=dip_image(ifftshift(double(in)));
else
    out=builtin('ifftshift',in);    % call the builtin version
end
