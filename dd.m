% res=dd(vargin)  : Display a series of matlab (or dip_image) datatype images as an N-dimensional dip_image object and returns
% varargin : series of images to join and display
% 
% Example
% a=readim;
% b=double(readim('orka'));
% dd a b
% joined=dd(a,b)

function res=dd(varargin)
maxdim=1;
for p=1:numel(varargin)
    if isa(varargin{p},'char')
        varargin{p}=evalin('caller',varargin{p}); % convert the string into the variable
    end
    if iscell(varargin{p})
        varargin{p}=dd(varargin{1}{:});
    end
    varargin{p}=dip_image(varargin{p});  % convert each argument to dip_image if necessary.
    sz=size(varargin{p});
    maxdim=max(maxdim,find(sz > 1,1,'last'));
end
fprintf('Displaying %d images concatinated along direction %d \n',numel(varargin),maxdim);
res=cat(maxdim+1,varargin{:});
