% [myPSF,FWHM,myResiduum,beadsAt,AllParam,subPixelShifts]=ExtractMultiPSF(myim,mysize,Pixelsize,doIndivFit,beadsAt,allowXY,subPixelShifts) :  Extracts mutliple PSFs (beads) and overlays them to obtain final PSF
% myim : input image to select the beads from
% mysize : size of the ROI to extract around each bead. This should be a vector
% Pixelsize : physical dimension of the pixel (can be a vector)
% doIndevFit  : fit each extracted PSF individually. All FWHMs and all the fit parameters are returned
% beadsAt : if not empty, the user is asked to click to select the bead positions
% allowXY: Flag. If 1, the fit will be independent in XY for FWHM: Default=0
% 
% Output:
% myPSF: The overlayed (summed) PSFs.
% FWHM:  the full width at half maximum as obtained from the Gaussian fit(s)
% myResiduum: The difference between fit and overlayed PSFs.
% beadsAt: The positions of all the clicks (precise to one pixel). Can be fed back into the algorithm.
% AllParam: A matrix with all the fit paramters
% subPixelShifts: The precise shifts used for the alignment between all the PSFs
%
% Example : [myPSF,FWHM,myResiduum,beadsAt,AllParam,subPixelShifts]=ExtractMultiPSF(myim,[16 16],[0.065 0.065],1,[])

function [myPSF,FWHMs,myResiduum,beadsAt,AllParam,subPixelShifts]=ExtractMultiPSF(myim,mysize,Pixelsize,doIndivFit,beadsAt,allowXY,subPixelShifts)
% beadsAt=[180 504; 207 489; 253 236; 399 276; 262 54; 364 49; 393 168; 410 156; 300 139; 190 138; 376 282; 182 459; 333 449; 367 467; 415 453; 204 550; 175 579];
if nargin < 7
    subPixelShifts={};
end
if nargin < 6
    allowXY=0;
end
if nargin < 4
    doIndivFit=0;
end
if nargin < 2
    mysize=[10 10];
end

if nargin < 3
    Pixelsize=95;
end

myres=myim;
if nargin < 5 || isempty(beadsAt)
    myres
    fprintf('Please select bead coordinates by left clicking. Finish with right click\n');
    beadsAt=dipgetcoords(100);
    fprintf('%d beads selected\n',size(beadsAt,1)-1);
    beadsAt(end,:)=[];
end

[sumBdsSum,mybdsSum,sumBdsSum2,mybdsSum2,subPixelShifts]=MeanFromCoord(myres,beadsAt,mysize,[],subPixelShifts);  % Extracts the beads and aligns them and calculates the sum image

%myresSum=para.res_sumIm;
%[sumBdsSum,mybdsSum]=MeanFromCoord(myresSum,beadsAt,s)

%im = readtimeseries('beads1_000');
%myres=squeeze(sum(im(:,:,3:5)));
%[sumBdsRaw,mybdsRaw]=MeanFromCoord(myres,beadsAt,s)

% [params,res,fitted]=FitDataNDFast([10 10 0; 1000 s/2 s/2],sumBdsSum,300,'idiv')
%if(ndims(sumBdsSum) > 2)   % This is 3D data
%    sumBdsSum2D=sumBdsSum(:,:,floor(size(sumBdsSum,3)/2));
%    mysize=mysize(1:2);
%else
%    sumBdsSum2D=sumBdsSum;
%end

%fprintf('Fitting PSF ... \n')
%[fittedSum,paramsSum,idiv,myfunct] = FitDataND('c(1)*exp(-((x{1}-c(2)).^2+(x{2}-c(3)).^2)/c(4))+c(5)',[max(sumBdsSum2D)-min(sumBdsSum2D) 0.1 0.1 (norm(mysize)/8)^2 min(sumBdsSum2D)],sumBdsSum2D,1000);
%FWHM=sqrt(paramsSum(4))*Pixelsize(1) * 2*sqrt(log(2))  % to get FWHM
if isempty(subPixelShifts)
    [FWHMs,myPSF,myResiduum,paramsSum] = FitNDPSF(sumBdsSum,Pixelsize,[],allowXY);
else
    [FWHMs,myPSF,myResiduum,paramsSum] = FitNDPSF(sumBdsSum,Pixelsize,[],allowXY,0); % No centering of PSF!
end
FWHMs=abs(FWHMs);
fprintf('Sum Beads, FWHM is ');
fprintf('%g ',FWHMs);
fprintf('\n');

fprintf('Sum Beads, Height is ');
fprintf('%g ',paramsSum(8));
fprintf('\n');

if doIndivFit
    numbds=size(mybdsSum,2);
    myResiduum=newim([size(myResiduum) numbds]);
    % paramsSum(1)=paramsSum(1)/size(mybdsSum,2);  % correct for the sum being brighter than the individual ones
    if ndims(sumBdsSum) == 2
        paramsSum(8)=paramsSum(8)/numbds;  % correct for the sum being brighter than the individual ones
        alpha=1/paramsSum(8);
    else
        paramsSum(11)=paramsSum(11)/numbds;  % correct for the sum being brighter than the individual ones
        alpha=1/paramsSum(11);
    end
    paramsSum(1)=alpha*paramsSum(1)/numbds;   % correct for the sum background being brighter than the individual ones    
    paramsSum(8)=1;  % because the bead will be individually multiplied by alpha
    paramsSave=paramsSum;
    for n=1:numbds
        myBd=alpha*mybdsSum{n};
        %if(size(myBd) > 2)   % This is 3D data
        %    myBd2D=myBd(:,:,floor(size(myBd,3)/2));
        %end
        
        fprintf('Fitting PSF number %d... \n',n)
        paramsSum=paramsSave;
        if isempty(subPixelShifts)
            [iFWHMs,imyPSF,imyResiduum,iParam] = FitNDPSF(myBd,Pixelsize,paramsSum,allowXY);
        else
            [iFWHMs,imyPSF,imyResiduum,iParam] = FitNDPSF(myBd,Pixelsize,paramsSum,allowXY,0);
        end
        % [fittedSumI,paramsSumI,idivI,myfunctI] = FitDataND('c(1)*exp(-((x{1}-c(2)).^2+(x{2}-c(3)).^2)/c(4))+c(5)',paramsSum,myBd,1000);
        AllParam{n}=iParam;AllParam{n}(end)=AllParam{n}(end)/alpha;AllParam{n}(1)=AllParam{n}(1)/alpha;
        iFWHMs=abs(iFWHMs);
        FWHM(n,:)=iFWHMs;  % to get FWHM
        if length(iFWHMs) > 2
            fprintf('Bd Number %d, FWHM is [%g,%g,%g]\n',n,iFWHMs);
        else
            fprintf('Bd Number %d, FWHM is %g,\n',n,iFWHMs);
        end
        if length(size(myResiduum)) == 2
            myResiduum(:,n-1)=imyResiduum;
        elseif length(size(myResiduum)) == 3
            myResiduum(:,:,n-1)=imyResiduum;
        elseif length(size(myResiduum)) == 4
            myResiduum(:,:,:,n-1)=imyResiduum;
        else
            error('Unsupported dimension number');
        end
    end
    if length(iFWHMs) > 2
        fprintf('FWHM of overlay is [%g %g %g], Mean if individual FWHMs is [%g %g %g], Stddev is [%g %g %g]\n',FWHMs,mean(FWHM),std(FWHM));
    else
        fprintf('FWHM of overlay is %g, Mean if individual FWHMs is %g, Stddev is %g\n',mean(FWHMs),mean(FWHM),std(FWHM));
    end
    FWHMs=FWHM;  % return the individual fit results
else
    AllParam=paramsSum;
end

%[fittedHRes,paramsHRes,idiv,myfunct] = FitDataND('c(1)*exp(-((x{1}-c(2)).^2+(x{2}-c(3)).^2)/c(4))+c(5)',[max(sumBdsHighRes) 0 0 20 min(sumBdsHighRes)],sumBdsHighRes,1000);
%FWHM=sqrt(paramsHRes(4))*Pixelsize* 2*sqrt(log(2))  % to get FWHM
%sumBdsHighRes-fittedHRes
