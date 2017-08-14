% riftshift : shifting the zero frequency back to the corner

function out = riftshift(in)
toshift=size(in)*0;
toshift(1)= - floor(size(in,1)/2);
out=circshift(in,toshift);