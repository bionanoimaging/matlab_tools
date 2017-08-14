% function out=funnyFunction(in) : A function 
% 
% Example:
% myWindow=funnyFunction(2.0*xx([256 1],'freq')) .* funnyFunction(2.0*yy([1 256],'freq'))

function out=funnyFunction(in)
out=in.*0;
in=abssqr(in);
mymask=in < 1.0;
out(mymask) = exp(-1/(1-in(mymask)));
