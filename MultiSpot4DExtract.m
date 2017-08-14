function [res4d, spotpos] = MultiSpot4DExtract5(img, correlroi, thresholdvalue, spotpos, basisvectors, zeroPos)
    % img: image to process
    % correlroi: already extracted ROI of a single spot which will be used for correlations or to determine the size of the sub ROI to extract around each spot
    % spotpos: if not empty ([]) it states the subpixel spot positions and the values. Will also be returned and can be used for successive calls. 
    % basisvectors: two vectors stating the approximate distance between the spots. 
    % thresholdvalue: A value to define the minimum quality of correlations to could as spot. 
    % res4d: 4-dimensional stack of extracted ROIs: (MillerindexX, MillerindexY, ROICoordinateX, ROICoordinateY)

if nargin<5
    %getting "img" (mydatbg)
    img1 = readim('2016-02-17 0001 1-30s - NIR LED 100mA + Raman Filter W.tiff'); %img1 = readim('2016-02-17\2016-02-17 0001 1-30s - NIR LED 100mA + Raman Filter W.tiff');
    img2 = readim('2016-02-17 0002 1-30s - NIR LED 100mA + Raman Filter W.tiff');
    img3 = readim('2016-02-17 0003 1-30s - NIR LED 100mA + Raman Filter W.tiff');
    img4 = readim('2016-02-17 0004 1-30s - NIR LED 100mA + Raman Filter W.tiff');
    mybg1 = readim('2016-02-17 BG01 1-30s - NIR LED 10mA W.tiff');
    mybg2 = readim('2016-02-17 BG02 1-30s - NIR LED 10mA W.tiff');
    mybg3 = readim('2016-02-17 BG03 1-30s - NIR LED 10mA W.tiff');
    mybg4 = readim('2016-02-17 BG04 1-30s - NIR LED 10mA W.tiff');
    img = (img1+img2+img3+img4)/4 - (mybg1+mybg2+mybg3+mybg4)/4;%mydatbg
    %getting "correlroi" (mytemp)
    correlroi = img(802:895,646:657);%data taken from image
    correlroi = extract(correlroi - min(correlroi)-300,size(img));
    %getting "thresholdvalue" (mythresh)
    thresholdvalue = 1e7;
    %getting "spotspos" (mymax)
    mycc = real(ift(ft(img) .* conj(ft(correlroi))));%Correlation
    mycc(mycc < thresholdvalue)=0;%no binary image
    [spotpos,spotposv]=findmaxima(mycc);%sub-pixel coordinates and correspondent interpolated values
        %From Subpixel to Pixel, building Spot Image
            SpotImage = newim(img);
            myidx  =floor(spotpos(:,1)) * size(mycc,2) + floor(spotpos(:,2)) ; %from index a(x,y) to index a(n)
            SpotImage(myidx)=spotposv;%SPOT IMAGE
    %getting basisvectors 
    zeroPos=[803 639];% taken from spot image
    basis1=([1073 714]-zeroPos)/6.0 * 1.1; % taken from spot image
    basis2=([747 920]-zeroPos)/5.0;% taken from spot image
    basisvectors = cat(3,basis1,basis2);%basis1 = basisvectors(:,:,1) and basis2 = basisvectors(:,:,2)
else
    mycc = real(ift(ft(img) .* conj(ft(correlroi))));%Correlation
    mycc(mycc < thresholdvalue)=0;%no binary image
    SpotImage = newim(img);
    myidx  =floor(spotpos(:,1)) * size(mycc,2) + floor(spotpos(:,2)) ; %from index a(x,y) to index a(n)
    SpotImage(myidx)=100;%SPOT IMAGE, just to observ spots
end


%Matrix problem: Miller Ccoordinates
M=[basisvectors(:,:,1);basisvectors(:,:,2)];%it was M=[basis1;basis2];
% MeasuredCoord = M * NewCoord
% NewCoord = M^-1 * MeasuredCoord
MillerIndexCoord = inv(M') * transpose(spotpos - repmat(zeroPos,[size(spotpos,1) 1]));
%MillerindexX = MillerIndexCoord(1,:);
%MillerindexY = MillerIndexCoord(2,:);

%ax, ay will be spot Images, instead of intensity values they have
%assigned MillerIndexCoorX and MillerIndexCoorY respectively
ax=newim(img);
ax(myidx)=MillerIndexCoord(1,:); % placing the MillerCoorX in its corresponding spot (in the spot Image)
ay=newim(img);
ay(myidx)=MillerIndexCoord(2,:); % placing the MillerCoorY in its corresponding spot (in the spot Image)

%Ranges of margings
myrngMinus=50; % this means plus minus
myrngPlus=100; % this means plus minus
 
mymask=SpotImage>0;%mask: binary image
%cleanning margings (left and right) space for placing spectral data (horizontal)
mymask(0:myrngMinus,:)=0;
mymask(size(mymask,1)-myrngPlus-1:end,:)=0;
ax(0:myrngMinus,:)=0;
ax(size(mymask,1)-myrngPlus-1:end,:)=0;
ay(0:myrngMinus,:)=0;
ay(size(mymask,1)-myrngPlus-1:end,:)=0;


% Ranges

    FullRngX = myrngMinus+myrngPlus+1;%length of spectra
    FullRngY = 10+1;

% Miller Coordinates   

    %mytemp = extract(repmat(dip_image([1:FullRng]),[1 1]),size(mycc))>0;
    mytemp = newim(1600,1200);
    mytemp(MidPosX(mytemp)-myrngMinus:MidPosX(mytemp)+myrngPlus,MidPosY(mytemp)-5:MidPosY(mytemp)+5)=1; % index for spectral coordinate %Index is not gradient
    %Assigning the values of each MillerCoor
    myXind = round(real(ift(ft(ax) .* ft(mytemp)))*sqrt(prod(size(mytemp))));% x Miller coordinates on each ax spot obtained by convolution of ax with spectral template
    myYind = round(real(ift(ft(ay) .* ft(mytemp)))*sqrt(prod(size(mytemp))));% y Miller coordinates  on each ay spot obtained by convolution of ax with spectral template
    
    fullMask = round(real(ift(ft(mymask) .* ft(mytemp)))*sqrt(prod(size(mytemp)))) > 0;% Mask of Spectra!
    
    % adjust the corner (adjusting intensities ONLY in regions of interest determined my fullMask)
    myXind(fullMask)=myXind(fullMask)- min(myXind(fullMask));
    myYind(fullMask)=myYind(fullMask)- min(myYind(fullMask));
    

%ROI Coordinates
    
    FullRngXIm = newim(FullRngX,1);%151 x 1
    FullRngXIm(0:end) = 0:FullRngX-1; %151 x 1
    FullRngXIm =  repmat(FullRngXIm,1,FullRngY);%151 x 11
   
    FullRngYIm = newim(1,FullRngY);% 1 x 11
    FullRngYIm(0:end) = 0:FullRngY-1; % 1 x 11
    FullRngYIm =  repmat(FullRngYIm,FullRngX,1);%1 x 11
    
    mytempX = newim(1600,1200);
    mytempX(MidPosX(mytempX)-myrngMinus:MidPosX(mytempX)+myrngPlus,MidPosY(mytempX)-5:MidPosY(mytempX)+5)=FullRngXIm; % index for spectral coordinate %index is a gradient
    mytempY = newim(1600,1200);
    mytempY(MidPosX(mytempY)-myrngMinus:MidPosX(mytempY)+myrngPlus,MidPosY(mytempY)-5:MidPosY(mytempY)+5)=FullRngYIm; % index for spectral coordinate %index is a gradient

    %ROICoordinateX
    mySindX = round(real(ift(ft(mymask) .* ft(mytempX)))*sqrt(prod(size(mytempX))));%convolution to have an index spectral coordinate on each maxima spot position
    %ROICoordinateY
    mySindY = round(real(ift(ft(mymask) .* ft(mytempY)))*sqrt(prod(size(mytempY))));
    %mySind: spectral index

%Preparing the 4D data matrix
   
    res4d=newim([max(myXind)+1,max(myYind)+1,max(mySindX)+1,max(mySindY)+1]);%4D
    myTotalInd4D = myYind + myXind*size(res4d,2) + mySindX*size(res4d,2)*size(res4d,1) + mySindY*size(res4d,3)*size(res4d,2)*size(res4d,1);
    cat(3,fullMask,myXind,myYind,mySindX,myTotalInd4D);
    
    %"Trick"
    tic
    %extracting the important information (from original image) which location is determined by binary fullMask
    res4d(double(myTotalInd4D(fullMask))) = img(fullMask);
    toc

 %Rescaling
    %allData2x=myrescaling2x_4D(allData4D)
    %ft4D=dip_fouriertransform(res4d, 'forward', [0 0 1 1]);
    %real(dip_fouriertransform(ft4D.*exp(1i*2*pi*ramp(ft4D,'freq')*0.5), 'inverse', [0 0 1 1]));
    
    
    
    
    