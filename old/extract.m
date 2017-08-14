% extract will extract (cut, pad) a specified region centered at a defined
% synopsis: out=extract(img,asize,center,value)
% coordinate and zero-pad if needed
% img: input image
% asize: size of output image
% center: center position (will be rounded)
% value: value to use for padding
% defaults:
%   center: floor(size(img)/2)
%   value: 0
%
%  Author: Rainer Heintzmann
function out=extract(img,asize,center,value)
isize=size(img);
if nargin < 4
    value=0;
end
if nargin < 3
    center=floor(isize/2);
else
end
idim=size(asize);
if idim(1) > 1
        asize = asize';
end
idim=size(center);
if idim(1) > 1
        center = center';
end
if length(asize) < length(isize)
    for d=length(asize)+1: length(isize)
        asize(d) = size(img,d);
    end
end
if length(center) < length(isize)
    for d=length(center)+1: length(isize)
        center(d) = floor(size(img,d)/2);
    end
end
N=ndims(img);
if isa(img,'cuda')
    global use_newim_cuda; tmp=use_newim_cuda; use_newim_cuda=1;
    out=newim_cuda(asize,datatype(img));   % problems with type conversion: +value;
    use_newim_cuda=tmp;
else
    out=newim(asize,datatype(img));   % problems with type conversion: +value;
end
if value~=0
    out=out+value;
end
srccenter=round(center');
srcstart=srccenter-floor(asize'/2);
srcend=srcstart+asize'-1;
dststart=zeros(N,1);
dststart(srcstart<0)=-srcstart(srcstart<0);
dstend=asize-1;
dstend(srcend>=isize')=dstend(srcend'>=isize)'-srcend(srcend'>=isize)+isize(srcend'>=isize)'-1;
srcend(srcend>=isize')=isize(srcend'>=isize)-1;
srcstart(srcstart<0)=0;
       
s = 'out(';
for ii=1:N
   s = [s 'dststart(' num2str(ii) '):dstend(' num2str(ii) '),'];
end
s(end) =[];
s = [s ')=img('];
for ii=1:N
   s = [s 'srcstart(' num2str(ii) '):srcend(' num2str(ii) '),'];
end
s(end) =[];
s = [s ');'];
eval(s); 
