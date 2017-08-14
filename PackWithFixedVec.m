% [newvec,numvars]=PackWithFixedVec(myvec,fixedvec,constvec,numdims) : creates a long vector with all entries using the description in fixedvec
%
% fixedvec defines how the variables should be treated: 
%   0 : variable in each of the Gaussian spots
%   1 : fixed value in each of the Gaussian spots
%   2 and higher : binds several instances of this same number together (global variable), this leads to a summing of the respective derivatives

function [newvec,numvars]=PackWithFixedVec(myvec,fixedvec,constvec,numdims)
skipdims=(2*numdims+1);
offset=1+numdims;  % backgournd + slopes
numSpots=(length(fixedvec)-offset)/skipdims;
if ~(numSpots == floor(numSpots))
    error('FixedMultiGaussDeriv: input vector has wrong length');
end

newvec=fixedvec';
newvec(fixedvec==0)=constvec;  % fill in the constants
numvars=sum(fixedvec==1);
newvec(fixedvec==1)=myvec(1:numvars);     % fill in the variables
for q=1:length(myvec)-numvars
    newvec(fixedvec== q+1)=myvec(numvars+q);     % fill in the globals
end

