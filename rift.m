% out=rift(in) : Simulates an inverse rft with a complex valued half complex array as input
function out=rift(in)
wasdip=isa(in,'dip_image');
in=double(rft2fft(in));
out=real(ifftn(in))*sqrt(prod(size(in)));
if wasdip
    out=dip_image(out);
end
