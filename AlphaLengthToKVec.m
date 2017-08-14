% KVec=GratingVecToKVec(Alpha,Length,ImgSize) : converts a an angle and distance of a real space grating to a KVector
% Alpha : Angle of the grating in real space in radiants. Can also be vectors, which will individually be processed
% Length : Grating perion. (= grating constant = distance from max to max). This can also be a vector, which is interpreted as different grating periods for each alpha.
% ImgSize : Size of the image (real space size = Fourier space size)
%
function KVec=AlphaLengthToKVec(Alpha,Length,ImgSize)
if numel(ImgSize) ~= 2
    error('AlphaLengthToKVec only works for 2-dimensional images (image size is different here)')
end
% ImgSize=floor(ImgSize/2)*2;  % adjusts for the behaviour of Fourier-transformations with uneven sizes.
% KVec = ([cos(Alpha') sin(Alpha')] .* repmat(ImgSize,[size(Alpha,2) 1])) ./ (repmat(Length.^2,[2 1]))';  % The transpose is needed to allow multiple vectors beeing processed at once

KVec = [cos(Alpha') sin(Alpha')] .* ((1./Length)' * ImgSize);
