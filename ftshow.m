% function img=ftshow(img,myexp)  : Fourier transforms data and returns the absolute magnitude of it to the power of myexp
% img : Image to Fourier transform
% myexp : exponent to use (default = 0.15)

function img=ftshow(img,myexp)
if nargin < 2
    myexp=0.15;
end
img =abs(ft(img)).^myexp;
end
