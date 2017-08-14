% out=rift2d(in) : performs a 2d inverse rft for each 2d slice in a dataset
function out=rift2d(in)
numdims=ndims(in);
if numdims > 4
    error('rift2d is only defined for up to 4 dimensions');
else
    in=expanddim(in,4);
end
out=newim(size(in,1),2*(size(in,2)-1),size(in,3),size(in,4),'single');
for e=0:size(in,4)-1
    for z=0:size(in,3)-1
        out(:,:,z,e)=rift(squeeze(in(:,:,z,e)));  % performs individual rift for each slice
    end
end
if numdims < 4
    mysize=size(out);
    out=reshape(out,mysize(1:numdims));
end
