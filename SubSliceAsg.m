% res=SubSliceAsg(img,myDim,myslice,val) : assign to a subslices in a dataset along a user-defined dimension
% img : dataset to subslice
% myDim : dimension along which to subslice
% myslice : one or multiple slicing position
%
% Example:
% a=readim('chromo3d');
% a=SubSliceAsg(a,3,[5],a(:,:,7))  % extract and display 3rd, 5th and 7th slice of the 3D data.
%
function res=SubSliceAsg(img,myDim,myslice,val)
if nargin < 3    
    myslice=MidPos(img,myDim);
end

    idx=num2cell(repmat(':',[1 ndims(img)]));
    idx{myDim}=myslice;
    S=struct('type','()','subs',{idx});
    res=subsasgn(img,S,val);
    