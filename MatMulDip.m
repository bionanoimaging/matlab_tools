% res=MatMulDip(myMat,mixImg,myDim) : Multiplies a matrix with a direction in a DipImage dataset
% myMat : matrix to multiply with
% myImg : image to multiply with
% myDim : dimension along which to apply the multiplication to. (default: last dimension = ndims)
%
% Example: unmixing
% MatMulDip([1 2;3 4],cat(3,readim,readim('orka')))
% 

function res=MatMulDip(myMat,myImg,myDim)
if nargin < 3
    myDim=ndims(myImg);
end

% This is the old code. IT DOES NOT WORK FOR ARRAY SIZES > 10 !!
% mixArray = newimar(size(myImg,myDim));
% for p=1:size(myImg,myDim)
%     mixArray{p}=(SubSlice(myImg,myDim,p-1));  % Iiiih: matlab based sub-slicing
% end
% unmixedArray = transpose(myMat) * mixArray;  % Pixelwise unmixing
% res=cat(myDim,unmixedArray{:});

% new code via matlab:
sz = size(myImg);
myImg=double(myImg);
szm = size(myImg);
szm(myDim)=size(myMat,1);
res = dip_image(reshape(transpose(myMat * transpose(reshape(myImg,[prod(sz(1:myDim-1)),sz(myDim)]))),szm));

