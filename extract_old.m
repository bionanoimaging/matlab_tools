% out=extract(img,asize,center,value) : will extract (cut, pad) a specified region centered at a defined. Also with padding/extrapolation.
% img: input image
% asize: size of output image (default:[] = original size)
% center: center position (will be rounded) (default : [] = center of image)
% value: value to use for padding (default: 0)
%        if 'cyclic' is specified for value a cyclic blending ist done
%        if 'gcyclic' is specified for value a cyclic blending ist done with Gaussian prefiltering
%
% Example:
% a  = readim('orka');
% q = extract(a,[300 300],[],'gcyclic')
% cat(1,ft(a),extract(ft(q),[256 256]))
%
%  Author: Rainer Heintzmann
function out=extract(img,asize,center,value)
global use_newim_cuda
isize=size(img);
if nargin < 4
    value=0;
end
if nargin < 3 || isempty(center)
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
    out=newim(asize,datatype(img));   % problems with type conversion: +value;
else
    tmp=use_newim_cuda;
    use_newim_cuda=0;
    out=newim(asize,datatype(img));   % problems with type conversion: +value;
    use_newim_cuda=tmp;    
end

srccenter=round(center');
if length(asize) > length(srccenter)
   srccenter(end+1:length(asize))=0;
   isize(end+1:length(asize))=1;
end

srcstart=srccenter-floor(asize'/2);
srcend=srcstart+asize'-1;
dststart=zeros(length(asize),1);
dststart(srcstart<0)=-srcstart(srcstart<0);
dstend=asize-1;
dstend(srcend>=isize')=dstend(srcend'>=isize)'-srcend(srcend'>=isize)+isize(srcend'>=isize)'-1;
srcend(srcend>=isize')=isize(srcend'>=isize)-1;
srcstart(srcstart<0)=0;

if isnumeric(value)  
    if value~=0
        out=out+value;  % the rest will be overwritten below
    end
end


s = 'out(';
for ii=1:ndims(out)
   s = [s 'dststart(' num2str(ii) '):dstend(' num2str(ii) '),'];
end
s(end) =[];
s = [s ')=img('];
for ii=1:ndims(img)
   s = [s 'srcstart(' num2str(ii) '):srcend(' num2str(ii) '),'];
end
s(end) =[];
s = [s ');'];
eval(s); 

if ~isnumeric(value)  % other case was dealt with above
    switch value
        case {'cyclic','gcyclic'}
            for d=1:ndims(out) % dimension
                pstart=dststart(d);
                pend=dstend(d);
                sges=pstart+(size(out,d)-pend);  % total size to fill
                eval(['p1=out' hyperplaneString(d,pstart,ndims(out)) ';']);
                eval(['p2=out' hyperplaneString(d,pend,ndims(out)) ';']);
                for p=0:pstart-1 % position
                    w=(p+(size(out,d)-pend))/sges;
                    w=(sin(w*pi/2))^2;
                    aline=p1*w+p2*(1-w);
                    if value(1)=='g'  % use gaussian blurring
                        aline=gaussf(aline,w*(1-w)/0.25);
                    end
                    eval(['out' hyperplaneString(d,p,ndims(out)) '=aline;']);
                    %eval(['out' hyperplaneString(d,p,ndims(out)) '=w;']);
                end
                for p=pend+1:size(out,d)-1 % position
                    w=(p-pend)/sges;
                    w=(sin(w*pi/2))^2;
                    aline=p1*w+p2*(1-w);
                    if value(1)=='g'  % use gaussian blurring
                        aline=gaussf(aline,w*(1-w)/0.25);
                    end
                    eval(['out' hyperplaneString(d,p,ndims(out)) '=aline;']);
                end
            end
        otherwise
            error('extract: Unknown blend method');
    end
        
end
    