% GratingVec=KVecToGratingVec(KVec,ImgSize) : Converts a vector in Fourier-space into a grating vector in real space.
% KVec : Vector in Fourier space (measured from the center zero-frequency pixel)
% GratingVec : Vector of the grating in real space.
% ImgSize : Size of the image (real space size = Fourier space size). Note: This routine will automatically account for uneven image sizes and their meaning in Fourier-space.

function GratingVec=KVecToGratingVec(KVec,ImgSize)
% ImgSize=floor(ImgSize/2)*2;  % adjusts for the behaviour of Fourier-transformations with uneven sizes.
myGratingVec_div_L2 = KVec ./ ImgSize;
GratingVec = myGratingVec_div_L2 ./ (myGratingVec_div_L2 * myGratingVec_div_L2');
