% [Alpha, Length]=KVecToAlphaLength(KVec,ImgSize,Pixelsizes) : Converts a vector in Fourier-space into a grating vector in real space.
% KVec : Vector in Fourier space (measured from the center zero-frequency pixel)
% GratingVec : Vector of the grating in real space.
% ImgSize : Size of the image (real space size = Fourier space size). Note: This routine will automatically account for uneven image sizes and their meaning in Fourier-space.
% Pixelsizes : optional argument with the pixelsizes, defining the units of the Length. default=[1 1 ..];

function [Alpha, Length]=KVecToAlphaLength(KVec,ImgSize,Pixelsizes)
GVec=KVecToGratingVec(KVec,ImgSize);
if nargin >= 3
    GVec=GVec.*Pixelsizes;
end
[Alpha,Length]=GratingVecToAlphaLength(GVec);

