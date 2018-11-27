% [spotIm,Meas,showIm]=DrawSpots(SpotList,ImSize,Bg) : draws Gaussian spots into an image. Poisson noise is also applied.
% SpotList: List of spot intensity, XY postions, widthX, widthY, .
% ImSize: Size of the output image to generate
% Bg: Background level to apply
% 
% Example:
% [spotIm,Meas,showIm]=DrawSpots(cat(2,rand(20,1)*100,rand(20,2)*1000,rand(20,2)*10),[1024 1024],10)
%
function [spotIm,Meas,showIm]=DrawSpots(SpotList,ImSize,Bg)
NumSpots = size(SpotList,1);
if numel(ImSize)==1
    ImSize=[ImSize ImSize];
end
spotIm = newim(ImSize)+Bg;
showIm = newim(ImSize)+Bg;
for j=1:NumSpots
    Intensity = SpotList(j,1);
    posX = SpotList(j,2);
    posY = SpotList(j,3);
    sigX = SpotList(j,4);
    sigY = SpotList(j,5);
    spot = Intensity*exp(- (((xx(ImSize,'corner') - posX)/sigX) .^2 + ((yy(ImSize,'corner') - posY)/sigY) .^2));
    showIm = showIm+Intensity*(((abs(xx(ImSize,'corner') - posX)<sigX) * (abs(yy(ImSize,'corner') - posY)<1.1))+(abs(xx(ImSize,'corner') - posX)<1.1) * (abs(yy(ImSize,'corner') - posY)<sigY));
    spotIm =spotIm + spot;
end
Meas = noise(spotIm,'poisson');

