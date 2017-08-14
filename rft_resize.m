% resized=rft_resize(myrft,afactor)  : expands or shrinks in Fourier-space. Useful for rft-based rescaling (like fft).
% myrft: input in Fourier-space
% afactor: factor to rescale by. This can be a scalar (applied to each dimension) or a vector with a factor for each dimension.
% 
% Example:
% a=readim('chromo3d');
% fa=rft(a)
% res=rift(rft_resize(fa,2))
% a2=rift(rft_resize(rft(res),0.5))

function resized=rft_resize(myrft,afactor)
if nargin < 2
    afactor=2;
end
if length(afactor) == 1
    afactor=repmat(afactor,[1 ndims(myrft)]);
end
oldsize=size(myrft);
midold=floor(oldsize/2);
newsize=floor(afactor.*oldsize);
if ndims(myrft)==1
    %if mod(oldsize(1),2)
        newsize(1)=floor(afactor(1).*(oldsize(1)-1))+1;
    %end
else
    %if mod(oldsize(2),2)
        newsize(2)=floor(afactor(2)*(oldsize(2)-1))+1;
    %end
end
if ndims(myrft)>3
    newsize(4)=oldsize(4);
end
resized=newim(newsize,'scomplex');
for d=1:length(oldsize)
    if oldsize(d)>newsize(d)
        oldsize(d)=newsize(d);
        midold(d)=floor(newsize(d)/2);
    end
end
midoldL=midold;
midoldL(mod(oldsize,2)==0)=midoldL(mod(oldsize,2)==0)+1;  % necessary, because even and odd size rfts are slightly different in the way they deal with the border pixels
    
if ndims(myrft)==1
    resized(0:oldsize-1)=myrft;  % rest remains zero
elseif ndims(myrft)==2
    resized(0:midoldL(1)-1,0:oldsize(2)-1)=myrft(0:midoldL(1)-1,0:oldsize(2)-1);  % left part
    resized(end-midold(1)+1:end,0:oldsize(2)-1)=myrft(end-midold(1)+1:end,0:oldsize(2)-1); % right part 
elseif ndims(myrft)==3
    resized(0:midoldL(1)-1,0:oldsize(2)-1,0:midoldL(3)-1)=myrft(0:midoldL(1)-1,0:oldsize(2)-1,0:midoldL(3)-1);  % left part
    resized(end-midold(1)+1:end,0:oldsize(2)-1,0:midoldL(3)-1)=myrft(end-midold(1)+1:end,0:oldsize(2)-1,0:midoldL(3)-1); % right part 
    resized(0:midoldL(1)-1,0:oldsize(2)-1,end-midold(3)+1:end)=myrft(0:midoldL(1)-1,0:oldsize(2)-1,end-midold(3)+1:end);  % left part
    resized(end-midold(1)+1:end,0:oldsize(2)-1,end-midold(3)+1:end)=myrft(end-midold(1)+1:end,0:oldsize(2)-1,end-midold(3)+1:end); % right part 
elseif ndims(myrft)==4  % only do it for each element
    for e=0:oldsize(4)-1
    resized(0:midoldL(1)-1,0:oldsize(2)-1,0:midoldL(3)-1)=myrft(0:midoldL(1)-1,0:oldsize(2)-1,0:midoldL(3)-1);  % left part
    resized(end-midold(1)+1:end,0:oldsize(2)-1,0:midoldL(3)-1)=myrft(end-midold(1)+1:end,0:oldsize(2)-1,0:midoldL(3)-1); % right part 
    resized(0:midoldL(1)-1,0:oldsize(2)-1,end-midold(3)+1:end)=myrft(0:midoldL(1)-1,0:oldsize(2)-1,end-midold(3)+1:end);  % left part
    resized(end-midold(1)+1:end,0:oldsize(2)-1,end-midold(3)+1:end)=myrft(end-midold(1)+1:end,0:oldsize(2)-1,end-midold(3)+1:end); % right part 
    end
else
    error('rft_resize: unsupported number of dimensions');
end
