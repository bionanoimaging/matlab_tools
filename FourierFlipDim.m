% res=FourierFLipDim(data,toflip) : flips multiple dimensions in Fourier-space and correct to maintain the center position independet of size
function data=FourierFLipDim(data,toflip)
mysize=size(data);
for d=1:length(toflip)
    if (toflip(d))
        if mod(mysize(d), 2) == 0  % even
            newmid=floor(mysize/2);
            newmid(d)=newmid(d)-1;            
            data=extract(flipdim(data,d),[],newmid);  % accounts for negative scalings, which are not covered by extract
        else
            data=flipdim(data,d);  % accounts for negative scalings, which are not covered by extract
        end
    end    
end
