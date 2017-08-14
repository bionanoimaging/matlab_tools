% Plots the contrast with respect to the radial position (= distance to the center of the image).
% Synopsis:
% [mycontr,myfreq]=RadialContrast(img,angfreq,color)
% img: image from which to plot the contrast
% angfreq: frequency of the sampling. Default: 21. This should correspond to the number of spokes divided by two
% color: color of the curve (syntax of the native Matlab plot function). Default: 'b' (blue)

function [mycontr,myfreq]=RadialContrast(img,angfreq,color)
if nargin < 2
    angfreq=21; % to correspond to the standard spokes object
end
if nargin < 3
    color='b'; % to correspond to the standard spokes object
end
angfreq=angfreq*2;

myrad=CartImg2Pol(img);

myfreq=dip_fouriertransform(myrad,'forward',[1 0]);

midf=floor(size(myfreq,1)/2);
contrast=2*abs(myfreq(midf+angfreq,:))/abs(myfreq(midf,:));

plot(0:size(contrast,2)-1,double(contrast),color);
title('Contrast Curve');
xlabel('Radius');
ylabel(sprintf('Contrast at %g cycles/ring',angfreq));

if nargout> 1
    mycontr=contrast;
end
