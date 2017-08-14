% [res,deriv]=FixedMultiGaussDeriv(myvec,fixedvec,numdims) : fixes a number of paramters according to user preferences and calls MultiGaussMSEDeriv
% fixedvec defines how the variables should be treated: 
%   0 : variable in each of the Gaussian spots
%   1 : fixed value in each of the Gaussian spots
%   2 and higher : binds several instances of this same number together (global variable), this leads to a summing of the respective derivatives

function [res,dervec]=FixedMultiGaussDeriv(myvec,fixedvec,constvec,numdims)
[newvec,numvars]=PackWithFixedVec(myvec,fixedvec,constvec,numdims);  % expands into a vector with all entries

[res,deriv]=MultiGaussMSEDeriv(newvec);

dervec=CondenseGradientWithFixedVec(deriv,fixedvec,numvars);
