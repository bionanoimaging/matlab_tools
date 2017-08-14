% [res,otfs]=multiConvolve(images,psfs,otfs) : Apply a convolution to a series of images stored in a cell array.
function [res,otfs]=multiConvolve(images,psfs,otfs)
numimages=numel(images);
numpsfs=numel(psfs);
for n=1:numpsfs
    if nargin < 3
        otfs{n}=ft(psfs{n});
    end
    res{n}=real(ift(ft(images{mod(n-1,numimages)+1}) .* otfs{n}));
end