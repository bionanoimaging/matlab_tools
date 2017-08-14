% rftshift : shifting the zero frequency of a rft routine to the center

function out = rftshift(in)
toshift=size(in)*0;
toshift(1)=floor(size(in,1)/2);
out=circshift(in,toshift);