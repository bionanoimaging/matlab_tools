function [spotIm,Meas,showIm]=DrawSpots(SpotList,ImSize,Bg)
NumSpots = size(SpotList,1);
spotIm = newim([ImSize ImSize])+Bg;
showIm = newim([ImSize ImSize])+Bg;
for j=1:NumSpots
    Intensity = SpotList(j,1);
    posX = SpotList(j,2);
    posY = SpotList(j,3);
    sigX = SpotList(j,4);
    sigY = SpotList(j,5);
    spot = Intensity*exp(- (((xx([ImSize ImSize],'corner') - posX)/sigX) .^2 + ((yy([ImSize ImSize],'corner') - posY)/sigY) .^2));
    showIm = showIm+Intensity*(((abs(xx([ImSize ImSize],'corner') - posX)<sigX) * (abs(yy([ImSize ImSize],'corner') - posY)<1.1))+(abs(xx([ImSize ImSize],'corner') - posX)<1.1) * (abs(yy([ImSize ImSize],'corner') - posY)<sigY));
    spotIm =spotIm + spot;
end
Meas = noise(spotIm,'poisson');

