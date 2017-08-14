% a=SumAbsSqrRFT(v) : calculates the sum of abssqr values in RFT Fourier space, such that it corresponds to sum(abssqr(ft(rift(v))). This is a bit more complicated.
function a=SumAbsSqrRFT(v)
%a=2*sum(abssqr(v(:,1:end)))+sum(abssqr(v(:,0)));
if ndims(v) == 2
    a=2*sum(abssqr(v(:,1:end-1)))+sum(abssqr(v(:,0)))+sum(abssqr(v(:,end)));
elseif ndims(v) == 3
    a=2*sum(abssqr(v(:,1:end-1,:)))+sum(abssqr(v(:,0,:)))+sum(abssqr(v(:,end,:)));
elseif ndims(v) == 1
    a=2*sum(abssqr(v(1:end-1)))+sum(abssqr(v(0)))+sum(abssqr(v(end)));
else
    error('only 2d and 3d SumAbsSqr for RFTs are implemented')
end