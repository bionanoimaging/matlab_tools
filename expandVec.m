% res=expandVec(aVec,NumDims,value) : determines the mid Z position of an image in the dipimage sense = floor(size/2)
% aVec: Vestor to append values to
% NumDims: To which size does it need to be expanded (it also cuts!)
% value: which values will be filled in (default=1.0)

function res=expandVec(aVec,NumDims,value)
if nargin < 3
    value=1;
end
res=[aVec ones(1,NumDims-length(aVec))*value];
res=res(1:NumDims);
