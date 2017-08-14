% KVec=GratingVecToKVec(GratingVec,ImgSize,Pixelsizes) : converts a grating vector to a vector in Fourier-space.
% GratingVec : Vector of the grating
% ImgSize : Size of the image (real space size = Fourier space size). Note: This routine will automatically account for uneven image sizes and their meaning in Fourier-space.
% Pixelsizes : optional argument with the pixelsizes, defining the units of the Length. default=[1 1 ..];

function KVec=GratingVecToKVec(GratingVec,ImgSize,Pixelsizes)
if nargin >= 3
    GratingVec=GratingVec ./ Pixelsizes;
end
% ImgSize=floor(ImgSize/2)*2;  % adjusts for the behaviour of Fourier-transformations with uneven sizes.
L2 = GratingVec * GratingVec';
KVec = GratingVec .* ImgSize / L2;
