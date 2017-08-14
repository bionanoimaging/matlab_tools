function [img, correlroi, spotpos, basisvectors, zeroPos] = SpotImage(img1, img2, img3, img4, mybg1, mybg2, mybg3, mybg4, thresholdvalue)
%Inputs
%img1, img2, img3 and img4: The spectra images
%mybg1, mybg2, mybg3 and mybg4: The backgroung images
%thresholdvalue: A value to define the minimum quality of correlation.
%zeroPos: origin of Miller Coordinate System

%Outputs
%img: background substracted image
%correlroi: image with already extracted ROI of a single spot which will be used for correlations or to determine the size of the sub ROI to extract around each spot
%spotpos: subpixel spot positions and the values. Will also be returned and can be used for successive calls. 
%basisvectors:two vectors stating the approximate distance between the spots. 

if nargin<9
    %getting "img" 
    img1 = readim('2016-02-17 0001 1-30s - NIR LED 100mA + Raman Filter W.tiff'); %img1 = readim('2016-02-17\2016-02-17 0001 1-30s - NIR LED 100mA + Raman Filter W.tiff');
    img2 = readim('2016-02-17 0002 1-30s - NIR LED 100mA + Raman Filter W.tiff');
    img3 = readim('2016-02-17 0003 1-30s - NIR LED 100mA + Raman Filter W.tiff');
    img4 = readim('2016-02-17 0004 1-30s - NIR LED 100mA + Raman Filter W.tiff');
    mybg1 = readim('2016-02-17 BG01 1-30s - NIR LED 10mA W.tiff');
    mybg2 = readim('2016-02-17 BG02 1-30s - NIR LED 10mA W.tiff');
    mybg3 = readim('2016-02-17 BG03 1-30s - NIR LED 10mA W.tiff');
    mybg4 = readim('2016-02-17 BG04 1-30s - NIR LED 10mA W.tiff');
    %getting "thresholdvalue" 
    thresholdvalue = 1e7;
end

img = (img1+img2+img3+img4)/4 - (mybg1+mybg2+mybg3+mybg4)/4;%mydatbg

%getting "correlroi" (mytemp)
%Select up and down corner of the correlroi
img
fprintf('Please click on the corresponding upper left corner of the ROI');
up_corner=dipgetcoords(1);
fprintf('Please click on the corresponding lower right corner of the ROI');
down_corner=dipgetcoords(1);
correlroi = img(up_corner(1):down_corner(1),up_corner(2):down_corner(2));
correlroi = extract(correlroi - min(correlroi)-300,size(img));

%getting "spotspos" (mymax)
mycc = real(ift(ft(img) .* conj(ft(correlroi))));%Correlation
mycc(mycc < thresholdvalue)=0;
    [spotpos,spotposv]=findmaxima(mycc);%sub-pixel coordinates and correspondent interpolated values
        %From Subpixel to Pixel, building Spot Image
            SpotImage = newim(img);
            myidx  =floor(spotpos(:,1)) * size(mycc,2) + floor(spotpos(:,2)) ; %from index a(x,y) to index a(n)
            SpotImage(myidx)=spotposv;%SPOT IMAGE

%getting basisvectors 
SpotImage
fprintf('Please click on the central spot of the image \n');
zeroPos=dipgetcoords(1);%Select the middle spot of the image (origing of Miller Coordinates)
fprintf('Please click on the spot of base vector 1 \n');
basis1=dipgetcoords(1);%Select the base vector 1
basis1=basis1-zeroPos; 
fprintf('Please click on the spot of base vector 2 \n');
basis2=dipgetcoords(1);%Select the base vector 2
basis2=basis2-zeroPos;
basisvectors = cat(3,basis1,basis2);%basis1 = basisvectors(:,:,1) and basis2 = basisvectors(:,:,2)


