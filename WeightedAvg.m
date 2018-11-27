% res=WeightedAvg(FTImgs,OTFs,mydim) : calculate the weighted averaging of a number of datasets given a number of OTFs.
% FTImgs: Input (Fourier transformed) images to calculate the weighted average from
% OTFs: Transfer functions (i.e. strength of each image)
% mydim: dimension along which the various images (and OTFs) are stacked (default is last dimension = ndims(Imgs)).
% res: Weighted average
%
% Example:
% 
%
function res=WeightedAvg(Imgs,OTFs,mydim)
if nargin < 3
    mydim=ndims(Imgs);
end

SumOTF = sum(OTFs,[],mydim);
SqrOTF = abssqr(OTFs);
SqrSumOTF = abssqr(SumOTF);

validmask = SqrOTF > 1e-5;

myfilter = SumOTF ./ OTFs; % SubSlice(OTF,mydim,0);
myfilter(~validmask)=0.0;  % 

weights = SqrOTF ./ SqrSumOTF;
weights = weights ./ (sum(weights,[],mydim) + eps);
weights(~validmask)=0.0;  % FIX THISSS!

res=sum(Imgs .* myfilter .* weights,[],mydim);

res=real(ift2d(res));
