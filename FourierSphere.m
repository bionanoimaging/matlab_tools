% FourierSphere(sz,d) : 3D Fourier transform of a solid sphere in 3D
% sz : size of Fourier space (3d)
% d : diameter of sphere in real space
% b(k) = 3h(pi k d)
% h(x) = sin(x)/x^3 - cos(x)/x^2
% Example:
% f=FourierSphere([100,100,100],80,[2.0,1.0,1.0]);
% ift(f)
function res=FourierSphere(sz,d, pixelsize)
if nargin > 2
    k = rrscale(sz,1.0./pixelsize./sz);
else
    k = rr(sz,'freq');
end
x = pi*k*d;
res=sin(x)/x^3 - cos(x)/x^2;
res(floor(sz(1)/2),floor(sz(2)/2),floor(sz(3)/2))=1/3.0;
res = 3 * res/ sqrt(prod(sz));  % makes the integral equal to one
