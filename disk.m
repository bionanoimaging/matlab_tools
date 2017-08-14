% res=disk(mysize,myoffset,myrad) : computes a circular disk at a given position and radius
% Example:
% disk([256 256],[100 10],4)
function res=disk(mysize,myoffset,myrad)
res=(xx(mysize)-myoffset(1)).^2+(yy(mysize)-myoffset(2)).^2;
res = res <= myrad*myrad;
