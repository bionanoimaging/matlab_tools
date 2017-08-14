function imgout=ft1d(imgin)  % Two dimensional Fourier transform

transvec=zeros(1,ndims(imgin));
transvec(1)=1;

imgout=dip_fouriertransform(imgin,'forward',transvec);
