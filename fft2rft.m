% out=fft2rft(in) : deletes half of the already ffted data to convert it into an rft
function out=fft2rft(in)

%safety unequal sizes
%if isa(in,'dip_image') && (mod(size(in,2),2) || size(in,2)==1) || ((~ isa(in,'dip_image')) && (mod(size(in,1),2) || size(in,1)==1))
if isa(in,'dip_image') && (mod(size(in,2),2) && ~size(in,2)==1) || ((~ isa(in,'dip_image')) && (mod(size(in,1),2) && ~size(in,1)==1))
     error('The fft2rft function only accepts even size along dim 1/2 (matlab/dipImage), as the rift of the result would yield a different size');
end

if ndims(in) == 2
    out=in(:,0:floor(size(in,2)/2));
elseif ndims(in) == 3
    out=in(:,0:floor(size(in,2)/2),:);
elseif ndims(in) == 1
    out=in(0:floor(size(in,1)/2));
else
    error('fft2rft only defined for 1d, 2d and 3d arrays.');
end
