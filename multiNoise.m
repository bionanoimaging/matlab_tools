% res=multiNoise(images,noisetype,maxval) : Applies noise (e.g. Poisson noise) to a series of images stored in a cell array. Normalisation is done from the maximum pixel of the whole series.
function res=multiNoise(images,noisetype,maxval)
fac=1;
numimages=numel(images);
if nargin == 3
    fac=0;
    for n=1:numimages
        fac=max(fac,max(images{n})/maxval);
    end
end
for n=1:numimages
    res{n}=noise(images{n}/fac,noisetype);
end