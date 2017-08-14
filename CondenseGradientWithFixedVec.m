% [newvec,numvars]=CondenseGradientWithFixedVec(deriv,fixedvec,numvars) : creates a condensed vector with all entries using the description in fixedvec
% This is used for the gradient and thus different entries for the same global fit are summed in the gradient (fixedvec=2 or higher)
%
% fixedvec defines how the variables should be treated: 
%   0 : variable in each of the Gaussian spots
%   1 : fixed value in each of the Gaussian spots
%   2 and higher : binds several instances of this same number together (global variable), this leads to a summing of the respective derivatives

function dervec=CondenseGradientWithFixedVec(deriv,fixedvec,numvars)
numglobals=max(fixedvec)-1;
dervec=zeros(numvars+numglobals,1);
dervec(1:numvars)=deriv(fixedvec==1);
for q=1:numglobals
    dervec(1:numvars)=deriv(fixedvec==1);
    dervec(numvars+q)=sum(deriv(fixedvec== q+1));     % sum all the partial derivatives (rule of computing d/dx f(x,x) = d/da f(a,b) + d/db f(a,b)
end

