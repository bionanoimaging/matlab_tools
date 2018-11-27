% [TimeInput,TimeCast,resCuda]=GenericDebugCheck(fkt,TypeString,varargin)  : checks if the performance of a function is identical between two datatypes
% fkt : function to test
% TypeString : string of the datatype to cast to
% varargin : matlab inputs which are converted to other datatype
% 
% Example:
% a= [1 2 3];b=[4 5 6];
% GenericDebugCheck('plus','dip_image',@double,dip_image(a),b);
% 
function [TimeInput,TimeCast,resVal]=GenericDebugCheck(fkt,TypeString,myCast,varargin)
argMatlab=cell(1,nargin-3);
Nout = nargout-2;
if Nout<= 0
    Nout=1;
end
inargs=varargin;

for n=1:numel(inargs)
    if (isa(inargs{n},TypeString))
        argMatlab{n}=myCast(inargs{n});  % convert to double or dip_image
    else
        argMatlab{n}=inargs{n};
    end
end
resMatlab=cell(1,Nout);
tic
eval(['[resMatlab{:}]=' fkt '(argMatlab{:});']);
TimeCast=toc;
%% Compare Matlab to CudaMat result
resType=cell(1,Nout);
% setCudaSynchronize(1);
tic
eval(['[resType{:}]=' fkt '(inargs{:});']);
TimeInput=toc;
% setCudaSynchronize(0); pause(1);

for n=1:Nout
    D=double(resMatlab{n}- myCast(resType{n}));
    D=norm(D(:));
    M1=max(abs(resMatlab{n}(:)));
    M2=max(abs(resType{n}(:)));
    M=max(M1,M2);
    fprintf('checked function %s type %s .. relative error is %d percent\n',fkt,TypeString,D/M*100);
    if (D/M > 0.01)  % one percent error leads to scream
        error('Function %s: detected disagreement of results!',fkt);
    end    
end

fprintf('checked function %s type %s .. relative error is %d percent\n',fkt,TypeString,D/M*100);


if Nout>0
    resVal=resType{:};
end
