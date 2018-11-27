% res=rro(asize,offset,varargin) : generates an N-dimensional distance map to the origin. Like the "rr" function but with offset
% asize : size of the output image
% offset : vector to the origin
% varargin : further arguments like 'corner' to be passed to the ramp function
function res=rro(asize,offset,varargin)
    asum=(ramp(asize,1,varargin{:})-offset(1)).^2;
    for d=2:length(offset)
        asum=asum+(ramp(asize,d,varargin{:})-offset(d)).^2;
    end
    res=sqrt(asum);
    