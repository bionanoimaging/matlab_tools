% function [FWHM]=FitNDPSF(img,pixelsize,startParams,allowXY)  : Fits a 2 or 3d Gaussian to the PSF image
% img : 2 or 3d image
% pixelsize: vector of pixel sizes (x,y,z)
% startParams: vector of starting parameters. Default: []
% allowXY: Flag. If 1, the fit will be independent in XY for FWHM
% allowShift: Flag. If 1, the PSF will also be centered, acconding to the fit result. Default: 1
% allowOffsetCor: Flag. If 1, the PSF will be offset-corrected using the fit values.
function [FWHMs,myPSF,myResiduum,paramsSum]=FitNDPSF(img,Pixelsize,startParams,allowXY,allowShift, allowOffsetCor)
if nargin < 5
    allowOffsetCor=1;
end
if nargin < 5
    allowShift=1;
end
mysize=size(img);
alpha=1.0;
if nargin < 4
    allowXY=0;
end
if nargin < 3  || isempty(startParams)
     gimg=gaussf(img);
     [maxv,maxpos]=max(gimg);
     [minv,minpos]=min(gimg);
     maxpos = maxpos -floor(mysize/2);
     alpha=1/(maxv-minv);  % This scaling constant has been introduced for numerical reasons. Otherwise for low values, the fit does not converge
     img=img*alpha;
    if ndims(img) == 2
        % startParams=[max(img)-min(img) 0.1 0.1 (norm(mysize)/8)^2 min(img)];
        startParams=[alpha*min(img) 0 0 maxpos(1) maxpos(2) (mysize(1)/2) (mysize(2)/2) 1];
    else
        % startParams=[max(img)-min(img) 0.1 0.1 0.1 (mysize(1)/8) (mysize(2)/8) (mysize(3)/8) min(img)];
        [mv,maxpos]=max(gaussf(img));
        maxpos = maxpos -floor(mysize/2);
        startParams=[alpha*min(img) 0 0 0 maxpos(1) maxpos(2) maxpos(3) (mysize(1)/4) (mysize(2)/4) (mysize(3)/4) 1];
    end
end
if ndims(img) == 2
    if allowXY
        fixedParams=[1        0 0  1   1         1             1          1];  % Determines what to fit
    else
        fixedParams=[1        0 0  1   1         2             2          1];  % Determines what to fit
    end
else
    if allowXY
        fixedParams=[1        0 0 0  1   1   1        1             1          1  1];
    else
        fixedParams=[1        0 0 0  1   1   1        2             2          1  1];
    end
end
if ndims(img) == 2
    fprintf('Fitting 2D PSF ... \n')
    if numel(Pixelsize) > 2
        fprintf('Warning: 3 pixelsizes are given but the image is only 2D. Ignoring 3rd entry\n');
        Pixelsize=Pixelsize(1:2);
    end
    % [fittedSum,paramsSum,idiv,myfunct] = FitDataND('c(1)*exp(-((x{1}-c(2)).^2+(x{2}-c(3)).^2)/c(4))+c(5)',startParams,img,1000);
    [paramsSum,idiv,fittedSum] = FitDataNDFast(startParams,img,ndims(img),1000,'mse',fixedParams);
    % FWHMs=sqrt(paramsSum(4))*Pixelsize(1) * 2*sqrt(log(2))  % to get FWHM
    FWHMs=abs(paramsSum(6:7)) .* Pixelsize * 2*sqrt(log(2));  % to get FWHM
    myshift=-paramsSum(4:5);
    if nargout > 1
        if(ndims(img) > 2)   % This is 3D data
            myshift(3)=0;
        end
        if allowShift
            myPSF=shift(img-paramsSum(1),myshift)/alpha;
        else
            myPSF=(img-paramsSum(1))/alpha;
        end
    end
    if nargout > 2
        myResiduum=(img-fittedSum)/alpha;
    end
    if nargout > 3
        paramsSum(8)=paramsSum(8)/alpha;
    end
else
    fprintf('Fitting 3D PSF ... \n')
    % [fittedSum,paramsSum,idiv,myfunct] = FitDataND('c(1)*exp(-(((x{1}-c(2))/c(5)).^2+((x{2}-c(3))/c(6)).^2+((x{3}-c(4))/c(7)).^2))+c(8)',startParams,img,1000);
    [paramsSum,idiv,fittedSum] = FitDataNDFast(startParams,img,ndims(img),1000,'mse',fixedParams);
    % FWHMs=paramsSum(5:7) .* Pixelsize * 2 * sqrt(log(2))  % to get FWHM
    FWHMs=abs(paramsSum(8:10)) .* Pixelsize * 2 * sqrt(log(2));  % to get FWHM
    myshift=-paramsSum(5:7);
    if allowOffsetCor
        myOffset = paramsSum(1);
    else
        myOffset = 0;
    end
    if nargout > 1
        if allowShift
            myPSF=shift(img-myOffset,myshift)/alpha;
        else
            myPSF=(img-myOffset)/alpha;
        end
    end
    if nargout > 2
        myResiduum=(img-fittedSum)/alpha;
    end
    if nargout > 3
        paramsSum(11)=paramsSum(8)/alpha;
    end
end
