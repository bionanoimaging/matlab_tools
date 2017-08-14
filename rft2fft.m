% out=rft2fft(in) : Converts an rft half back to a full sized fft
function out=rft2fft(in)
if (~isa(in,'dipimage'))
    in=dip_image(in);
end
if ndims(in) == 2
    flipped=conj(flipdim(flipdim(in(:,1:end-1),2),1));
    flipped=cat(1,flipped(end,:),flipped(0:end-1,:));
out=cat(2,in,flipped);
elseif ndims(in) == 3
    flipped=conj(flipdim(flipdim(flipdim(in(:,1:end-1,:),2),1),3));
    flipped=cat(1,flipped(end,:,:),flipped(0:end-1,:,:));
    if (size(in,3)>1)
        flipped=cat(3,flipped(:,:,end),flipped(:,:,0:end-1));
    end
out=cat(2,in,flipped);
elseif ndims(in) == 1
    flipped=conj(flipdim(in(1:end-1),1));
    flipped=cat(1,flipped(end),flipped(0:end-1));
    out=cat(1,in,flipped);
else
    error('rft2fft only defined for 1d, 2d and 3d arrays.');
end
