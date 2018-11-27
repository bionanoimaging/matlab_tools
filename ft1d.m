% imgout=ft1d(imgin,mydir) : One dimensional Fourier transform,
% imgin : image to transform
% mydir: (default=1) direction to transform
%
function imgout=ft1d(imgin,mydir)  % One dimensional Fourier transform
if nargin < 2
    mydir=1;
end

transvec=zeros(1,ndims(imgin));
transvec(mydir)=1;

imgout=dip_fouriertransform(imgin,'forward',transvec);
