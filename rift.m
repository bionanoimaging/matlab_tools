% out=rift(in) : Simulates an inverse rft with a complex valued half complex array as input
function out=rift(in)
wasdip=isa(in,'dip_image');
in=double(rft2fft(in));
if wasdip
    out=real(ifftn(in))*sqrt(prod(size(in)));
    out=dip_image(out);
else
    out=real(ifftn(in));
end
