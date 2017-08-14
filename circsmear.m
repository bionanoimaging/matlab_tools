% res=circsmear(mycurve,asize)   creates a 2d image from a curve. E.g. to convert an OTF curve into a psf
% See also: radialmean
function res=circsmear(mycurve,asize)

myrr = rr(asize);
rrb = floor(myrr);
rrw = (myrr-rrb);
imb=msr2obj(dip_image(rrb+1,'uint16'),double(mycurve));
imt=msr2obj(dip_image(rrb+2,'uint16'),double(mycurve));
res=imb*(1-rrw)+imt*rrw;
