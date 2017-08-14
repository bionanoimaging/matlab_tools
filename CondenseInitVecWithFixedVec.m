% [condensed,constvec]=CondenseInitVecWithFixedVec(initvec,fixedvec) : creates a condensed vector with all entries using the description in fixedvec
% This is used for the initialisation vector. The first entry in the given full vector defines the value to be used
% initvec : initialisation vector giving one value to each variable (even of condensed)
%
% fixedvec defines how the variables should be treated: 
%   0 : variable in each of the Gaussian spots
%   1 : fixed value in each of the Gaussian spots
%   2 and higher : binds several instances of this same number together (global variable), this leads to a summing of the respective derivatives
%
% condensed : initialisation vector appropriately condensed (to be used for minfunc)
% constvec : vector to be bound to the calling function

function [condensed,constvec]=CondenseInitVecWithFixedVec(initvec,fixedvec)
constvec=initvec(fixedvec==0);  % collect all the constants
numvars=sum(fixedvec==1);
numglobals=max(fixedvec)-1;
condensed=zeros(1,numvars+numglobals);
condensed(1:numvars)=initvec(fixedvec==1);  % write the local variables first
for q=1:numglobals
    p=find(fixedvec== q+1);
    condensed(numvars+q)=initvec(p(1));     % the first appearance defines the constant or 
end
