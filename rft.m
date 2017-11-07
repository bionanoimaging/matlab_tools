% out=rft(in) : Simulates an rft with a dipimage (or matlab array) as an input object
function out=rft(in)

%safety unequal sizes
% if isa(in,'dip_image') && (mod(size(in,2),2) || size(in,2)==1) || ((~ isa(in,'dip_image')) && (mod(size(in,1),2) || size(in,1)==1))
if isa(in,'dip_image') && (mod(size(in,2),2) && ~size(in,2)==1) || ((~ isa(in,'dip_image')) && (mod(size(in,1),2) && ~size(in,1)==1))
     error('The rft function only accepts even size along dim 1/2 (matlab/dipImage), as the rift of the result would yield a different size');
end
if ndims(in) > 3
    error('rft only defined for 2d and 3d arrays.');
end

% in=readim('orka'); %test object (display the ft in log stretch)
b=fftn(double(in)); % dip_image Fourier Transform: whole plane

if isa(in,'dip_image') 
    b=dip_image(b)/sqrt(prod(size(b)));
else
    b=dip_image(b);
end
out=fft2rft(b);
if ~isa(in,'dip_image')
    out=double(out); % convert back to double
end
