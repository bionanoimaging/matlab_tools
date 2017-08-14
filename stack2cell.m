% res=stack2cell(mydat,mydim,dosqueeze) : extracts subslices from a stack and assembles them into a cell array
% mydat : data to extract slices from
% mydim : dimension at which slices are extracted (default: last dimension)
% dosqueeze : perform a squezze on the extracted subslices (default: 0)
%
% Exampe:
% a=readim;
% b=stack2cell(a,2);
% c=cat(2,b{:})
function res=stack2cell(mydat,mydim,dosqueeze)
mdim=ndims(mydat);
if nargin < 3
    dosqueeze=0;
end
if nargin < 2 || isempty(mydim)
    mydim=mdim;
end
if mydim > mdim
    error('stack2cell: Slice dimension to extract is bigger than number of available dimensions')
end
exsize=size(mydat,mydim);
margs=cell(1,mdim);
res=cell(1,exsize);
for p=0:exsize-1
    for d=1:mdim
        if (d==mydim)
            margs{d}=p;
        else
            margs{d}=':';
        end
    end
    if dosqueeze
        res{p+1}=squeeze(mydat(margs{:}));
    else
        res{p+1}=mydat(margs{:});
    end
end
