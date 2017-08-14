ScaleX=43;
h=kSimPSF({'na',1.49;'lambdaEm',645;'ri',1.52;'sZ',1;'scaleX',ScaleX;'scaleY',ScaleX});
im=readim * 0;
im(floor(rand(10,1)*prod(size(im))))=1;

imm=noise(1000*convolve(im,h),'poisson');
[myPSF,FWHM]=ExtractMultiPSF(imm,40,ScaleX)

imm=read3dtiff('reseau3-cleanedreduit.tif')
i1=squeeze(imm(:,:,0));
[myPSF,FWHM]=ExtractMultiPSF(i1,40,ScaleX)
cPSF=myPSF-mean(myPSF(end,:));
cPSF(cPSF<0)=0

ic=kcorrelator(DampEdge(imm(:,:,1:30),0.2,2,1,2),{'gf',0.1;'p',1;'plane',1;'fixplane',1});
icm=squeeze(mean(ic,[],3));
[myPSF,FWHM]=ExtractMultiPSF(icm,40,ScaleX)
%cPSF=myPSF-mean(myPSF(end,:));
cPSF=myPSF;
cPSF(cPSF<0)=0

tiffwrite('MyPSF.tif',cPSF,'yes')

