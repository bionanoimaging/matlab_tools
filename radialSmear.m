% function radialSmear(PSF0,XYSize,overSampling): smears an RZ series of radial curves into a stack of XY images by using linear interpolation
% RZImage: Series of radial curves, oversampled by overSampling along the radial direction
% XYSize: Size of the result image. The datatype will be equal to the input datatype of the RZImage.
% overSampling: Oversampling factor
%
% Example:
% a=rr()
% a=a(128:end,:)
% res=radialSmear(a+0.0,[200,256],1.0)
%
% Authors: Rainer Heintzmann, Dina Miora Ratsimandresy
%
function res=radialSmear(RZImage,XYSize,overSampling)
RZImage = permute(RZImage,[2 1]);
RZImage = double(RZImage);
tic
[SizeR,SizeZ] = size(RZImage);
im = double(rr(XYSize)/overSampling);
XYSize = [XYSize(2),XYSize(1)];
% below are the interpolators each as XY distances
interpXYLow = floor(im);
Idx1DLow = interpXYLow(:);
maxIdx = SizeR-2;
interpXYWeight = im(:) - Idx1DLow;
WeightsHigh = interpXYWeight(:);
Overflowmask = Idx1DLow > maxIdx;
WeightsHigh(Overflowmask) = 1.0;
WeightsLow = 1.0 - WeightsHigh;
Idx1DLow(Overflowmask) = maxIdx;
Idx1DLow = uint64(Idx1DLow);
% res = newim([XYSize,SizeZ],datatype(RZImage));
toc
res=cell(1,SizeZ);
for z=1:SizeZ
%    DstOffset = prod(XYSize)*z;
    SrcOffset = SizeR*(z-1)+1;
    CurIdx = SrcOffset+Idx1DLow;
    thisCurve = RZImage(CurIdx).*WeightsLow+RZImage(CurIdx+1).*WeightsHigh;
    % res(prod(XYSize)*z:prod(XYSize)*(z+1)-1) = thisCurve;
    res{z} = thisCurve;
end
toc
res = dip_image(reshape(cat(1,res{:}), [XYSize,SizeZ]));
toc
