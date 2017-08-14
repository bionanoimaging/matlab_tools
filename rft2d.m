% out=rft2d(in) : performs a 2d rft for each 2d slice in a dataset
function out=rft2d(in)
numdims=ndims(in);
if numdims > 4
    error('rift2d is only defined for up to 4 dimensions');
else
    in=expanddim(in,4);
end
out=newim(size(in,1),floor(size(in,2)/2)+1,size(in,3),size(in,4),'scomplex');
for e=0:size(in,4)-1
    for z=0:size(in,3)-1
        out(:,:,z,e)=rft(squeeze(in(:,:,z,e)));  % performs individual rift for each slice
    end
end

if numdims < 4
    mysize=size(out);
    out=reshape(out,mysize(1:numdims));
end
