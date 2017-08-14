% res=MatMulDip(myMat,mixImg,myDim) : Multiplies a matrix with a direction in a DipImage dataset
% myMat : matrix to multiply with
% myImg : image to multiply with
% myDim : dimension along which to apply the multiplication to. (default: last dimension = ndims)
%
% Example: unmixing
% 
% 

function res=MatMulDip(myMat,myImg,myDim)
if nargin < 3
    myDim=ndims(myImg);
end
mixArray = newimar(size(myImg,myDim));
for p=1:size(myImg,myDim)
    mixArray{p}=(SubSlice(myImg,myDim,p-1));  % Iiiih: matlab based sub-slicing
end

unmixedArray = transpose(myMat) * mixArray;  % Pixelwise unmixing

res=cat(myDim,unmixedArray{:});
