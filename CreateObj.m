function obj=CreateObj(mysize,numspikes)
obj=newim(mysize);
obj(rr(mysize,'freq')<0.4)=1;
obj(mod(phiphi(mysize)+pi,2*pi/numspikes)<2*pi/numspikes/2)=0;
