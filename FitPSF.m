% [FWHM]=FitPSF(img) : Fits a PSF and returns the Full width at half maximum value. Also usefule for chromatic shift maps
% Example:
% a=readim
% t=a>140
% extrapolate(a,t)
function [FWHM]=FitPSF(img, scales)
%[fittedSum,paramsSum,idiv,myfunct] = FitDataND('c(1)*exp(-((x{1}-c(2)).^2+(x{2}-c(3)).^2)/c(4))+c(5)',[(max(img)-min(img)) 1 1 8 min(img)],img,1000);
%FWHM=sqrt(paramsSum(4))* 2*sqrt(log(2));  % to get FWHM
% [paramsSum,mse,fittedSum] = FitDataNDFast([min(img) 5 (max(img)-min(img)) 0 0],img,2,1000,'mse');
if nargin < 2
    scales=[1 1];
end
img=squeeze(img)/max(img);
[paramsSum,mse,fittedSum] = FitDataNDFast([min(img) 0 0 0 0 3 3 (max(img)-min(img))],img,2,1000,'mse');
FWHM=paramsSum(6:7)* 2*sqrt(log(2)) .* scales;  % to get FWHM

img-fittedSum

xpos=(size(fittedSum,1)/2)+round(paramsSum(4));
ypos=(size(fittedSum,2)/2)+round(paramsSum(5));
%ypos=floor(size(fittedSum,2)/2)+round(paramsSum(3));
xmin=floor(xpos-5*FWHM(1)./scales(1));
xmax=ceil(xpos+5*FWHM(1)./scales(1));
ymin=floor(ypos-5*FWHM(2)./scales(2));
ymax=ceil(ypos+5*FWHM(2)./scales(2));

xtics=(([xmin:xmax]-xpos)*scales(1));
ytics=(([ymin:ymax]-ypos)*scales(2));
figure;

plot(xtics,fittedSum(xmin:xmax,round(ypos)),'b');hold on; plot(xtics,img(xmin:xmax,round(ypos)),'r');
plot(ytics,fittedSum(round(xpos),ymin:ymax),'c');plot(ytics,img(round(xpos),ymin:ymax),'m');
%dipshow(11,img-fittedSum); drawnow;
xlabel('position [µm]');
ylabel('intensity [au]');
legend({'x-cut','y-cut'});
