% [res,resFT,validmask]=WeightedAvgOTF(OTFs,FTImgs,mydim,myCOV,killHF) : calculate the weighted averaging of a number of datasets given a number of OTFs.
% OTFs: individual Transfer functions (assuming equal noise in the corresponding images)
% FTImgs: if not empty [] Fourier-transformed images will be subjected to weigthed averaging
%         if empty, only the noise-normalized OTF is computed instead
% mydim: dimension along which the various images (and OTFs) are stacked (default is last dimension = ndims(Imgs)).
% myCOV: If empty, the zero frequency strengths of the OTFs will each be used for noise normalization. If not empty, the OTFs are assumed als already scaled appropriately. If the covariance matrix is given here, it is used for an optimized weighted averaging.
%        The covariance is always assumed to be frequency independent, but can have a convariances between various OTFs.
%        If you do not want automatic scaling (e.g. in SIM), set myCOV to 1.0.
% killHF: If selected (1) the part beyond the passband will be set to zero. If killHF==0 the zero frequency weights will be used for this region
% 
% res: Weighted average image (or PSF_nn), normalized for flat noise
% resFT: Weighted average Fourier transform of image (or OTF_nn), normalized for flat noise
% validmask: 
%
% As an input this routine expects OTFs as obtained by ordinary image detection systems. This means
% that the noise level is assumed to be uniform over spatial frequency.
% 
% Example:
% 
%
function [res,resFT,validmask]=WeightedAvgOTF(OTFs,FTImgs,mydim,myCOV,killHF,ValidLimit)
if nargin < 6
    ValidLimit=1e-7;
end

if nargin < 5
    killHF=0;
end
if nargin < 3 || isempty(mydim)
    mydim=ndims(OTFs);
end
if nargin < 4
    myCOV=[];
end

% midX=MidPosX(OTFs); midY=MidPosY(OTFs);
if ~isempty(myCOV) 
    SigmaSqr = 1.0;
else
    SigmaSqr = abs(SubSlice(OTFs,mydim)); % takes the middle by default
    OTFs = OTFs ./ sqrt(SigmaSqr);   % Normalize all OTFs to have uniform variance. E.g. a second OTF which has a small value at zero transfers only few photons and therefore needs to be scaled up to weigh more.
end
if nargin < 2 || isempty(FTImgs)
    FTImgs=OTFs;
else
    if isempty(myCOV) || (prod(size(myCOV)) ==1)
        FTImgs = FTImgs ./ sqrt(SigmaSqr);   % Normalize also the images to have uniform variance (except for a general factor)
    else
        error('not implemented fully')
    end
end
if isreal(FTImgs)
    error('WeightedAvgOTF needs the Fourier-transform of the images not the images.')
end
if isreal(OTFs)
    error('WeightedAvgOTF needs the OTFs not the PSFs.')
end
%% at this point the OTF and/or image should have been normalized to uniform variance

% calculate the weights
if isempty(myCOV) || (prod(size(myCOV)) ==1)
    weights = conj(OTFs);
else
    invCOV = inv(myCOV); 
    COVWeights = sum(invCOV) ./ (sum(sum(invCOV)));
    myShape=size(OTFs).*0+1; myShape(mydim)=size(OTFs,mydim);
    weights = conj(OTFs) .* reshape(dip_image(COVWeights),myShape);
end

SumWgt = (sum(abssqr(OTFs),[],mydim) + eps);
weights = weights ./ sqrt(SumWgt);                % implements the final equation: sum(Img h* / sqrt(sum(abssqr(h)))    

validmask = SumWgt > (ValidLimit * max(SumWgt));  % abs(OTFs) > (1e-3 * max(abs(OTFs))) & 
if (killHF)
    notvalidmask = repmat(~validmask,[1 1 1 size(weights,4)]);
    weights(notvalidmask)=0.0; 
else
    notvalidmask=~validmask;
    for w=0:size(weights,mydim)-1
        wp=SubSlice(weights,mydim,w);  % (:,:,:,w);
        wp(notvalidmask)=MidVal(wp); % replaces the invalid stuff with the weight at the center. % wp(midX,midY,midZ,:); 
        weights=SubSliceAsg(weights,mydim,w,wp);
    end
end

resFT=squeeze(sum(FTImgs .* weights,[],mydim));  % Filter the OTF by the calculated filter

res=real(ift(resFT));
