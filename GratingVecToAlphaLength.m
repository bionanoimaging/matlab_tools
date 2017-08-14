% [Alpha Length] =GratingVecToAlphaLength(GratingVec) : converts a an angle and distance of a real space grating to a KVector
% GratingVec : Vector of the grating
% Alpha : Angle of the grating in real space in radiants
% Length : Grating vector length (= grating constant = distance from max to max)

function [Alpha,Length] =GratingVecToAlphaLength(GratingVec)
if numel(GratingVec) ~= 2
    error('GratingVecToKVec only works for 2-dimensional vectors')
end
Alpha = atan2(GratingVec(2),GratingVec(1));
Length=norm(GratingVec);
