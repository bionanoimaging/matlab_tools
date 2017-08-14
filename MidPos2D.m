% res=MidPos2D(anImg)  : determines the mid position of an image in the dipimage sense = floor(size()/2)

function res=MidPos2D(anImg)
res=floor(size(anImg)/2);
res=res(1:2);
