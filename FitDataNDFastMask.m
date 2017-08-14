% FitDataNDFastMask : Fits several ND Gaussians to N-dimensional image data
% Synopsis: [fitted, params, myfunct] = FitDataNDFastMask(InitParm, ParmMask, Data, numdims, maxiter, method, centerzero, ROISize)
% This function is superseeded by FitDataNDFast
%
%
% AFunctString should use the variable 'x' as a vector and x{1}, x{2} to 
% access its components, and the values to fit: c(1),c(2),...
%
% InitParm : Vector of initial parameters first row is global parameters,
% meaning is as follows: [global background, global width, intensity, posy, posx, intensity2, posy2, posx2...]
% other rows are local parameters
% ParmMask: is the mask deciding which parameters are fixed and which are
% variable in the optimization. It is a vector of same size as InitParm.
% (0) stands for variable parameter, 
% (1) stands for fixed parameter, and
% (2) stands for a global variable parametr. The first value in the ParmMask
% decides whether the background offset is fixed (1) or variable (0 or 2).
% Data : Experimental multidimensional (e.g. image) data to be fitted
% maxiter : maximum number of iterations in fit (default = 300)
% method : figure of merit to use for fitting 'mse', 'idiv' or 'fidiv'. 'mse' is least squares (assuming Gaussian noise),'idiv' is Czesar's i-divergence (Sterling approximation) and 'fidiv' is a fast version omitting data-dependent constants
% centerzero : zero will be corresponding to the center of the coordinate system, default: 0 = no zero centering
% ROISize : size of the ROI to use
%
% Example:
% a=noise(51.2*exp(-((xx(20,20,'corner')-8).^2/(3^2)+(yy(20,20,'corner')-11.2).^2/(4^2)))+33,'poisson');
% b=noise(22.2*exp(-((xx(20,20,'corner')-14.5).^2/(3^2)+(yy(20,20,'corner')-6.2).^2/(4^2))),'poisson');
% a+b
% [params, res, fitted, residual] = FitDataNDFastMask([10 30 8 10 2 2], [0 0 0 0 0 0],a+b, 2, 300, 'idiv')
% [params, res, fitted, residual] = FitDataNDFastMask([10 30 8 10 2 2 30 12 8 2 2], [0 0 0 0 2 2],a+b, 2, 300, 'idiv')
% [params, res, fitted, residual] = FitDataNDFastMask([30 30 8 10 2 2 30 12 8 2 2], [0 0 0 0 2 2],a+b, 2, 300, 'idiv',0,[8 8])
% fitted is the residual image (depending on method)
% params contains the result of the fit

function [params, res, fitted, residual] = FitDataNDFastMask(InitParm, ParmMask, Data, numdims, maxiter, mymethod, centerzero, ROISize)
global ParmVarIdxFull;
global ParmConstIdxFull;
global ParmGlobIdx;
global ParmGlobIdxFull;
global ParmConsts;
global nrspots;
global GlobVarBack;

if nargin < 4
    numdims = 2;
end
if nargin < 5
    maxiter = 3000;
end

if nargin < 6
    mymethod = 'mse';
end

if nargin < 7
    centerzero = 0;
end

if nargin < 8
    ROISize = [0 0 0];
end

if length(ParmMask) ~= 2+2*numdims % length(ParmMask)
    error('Parameter Mask needs to have the length 6 [GlobalBG, Intenensity, PosX, PosY, ..., SigmaX, SigmaY, ...]');
end
%if ~islogical(ParmMask)
%    ParmMask = logical(ParmMask);
%end

% Decide if the background should be a global variable or fixed
GlobVarBack = 0;
if ParmMask(1) == 0 || ParmMask(1) == 2
    tovary = InitParm(1);
    ParmConst = [];
    GlobVarBack = 1;
elseif ParmMask(1) == 1
    ParmConst = InitParm(1);
    tovary = [];
end
InitParm(1) = [];
ParmMask(1) = [];

% get indices of Variable, Global, and Constant Parameters
ParmVarIdx = find(ParmMask == 0);   % This is a list of indices pointing to the parameters which are supposed to be varied
ParmConstIdx = find(ParmMask == 1); % This is a list of indices pointing to the parameters which are constant
ParmGlobIdx = find(ParmMask == 2);  % This is a list of indices pointing to the parameters whcih are global variables (i.e. linked together and managed by only one variable parameter)

% get number of Gaussians to be fitted
nrspots = length(InitParm) / (2 * numdims + 1);
if nrspots*(2*numdims+1) ~= length(InitParm)
    error('Wrong number of initialization parameters. Needs to be global background and [I,px,py,sigmax,sigmay] for each spot, or the equivalent in 3D.');
end
% fprintf('Nr of spots to fit: %d\n',nrspots);

TV = zeros(1, length(ParmGlobIdx));
for i = 1 : length(ParmGlobIdx)
    TV(i) = mean(InitParm((2 * numdims + 1) * (0 : (nrspots - 1)) + ParmGlobIdx(i)));
end
tovary = horzcat(tovary, TV);

% Get the index of the global variables in the first spot
if ~isempty(ParmGlobIdx)
    ParmGlobIdxFull = repmat(ParmGlobIdx, nrspots, 1) + repmat((2 * numdims + 1) * (0 : (nrspots - 1))', 1, length(ParmGlobIdx));
    ParmGlobIdxFull = sort(ParmGlobIdxFull(:));
else
    ParmGlobIdxFull = ParmGlobIdx;
end
if ~isempty(ParmVarIdx)
    ParmVarIdxFull = repmat(ParmVarIdx, nrspots, 1) + repmat((2 * numdims + 1) * (0 : (nrspots - 1))', 1, length(ParmVarIdx));
    ParmVarIdxFull = sort(ParmVarIdxFull(:));
else
    ParmVarIdxFull = ParmVarIdx;
end
if ~isempty(ParmConstIdx)
    ParmConstIdxFull = repmat(ParmConstIdx, nrspots, 1) + repmat((2 * numdims + 1) * (0 : (nrspots - 1))', 1, length(ParmConstIdx));
    ParmConstIdxFull = sort(ParmConstIdxFull(:));
else
    ParmConstIdxFull = ParmConstIdx;
end
ParmConsts = horzcat(ParmConst, InitParm(ParmConstIdxFull));

tovary = horzcat(tovary, InitParm(ParmVarIdxFull))';


%if size(ParmMask) ~= size(InitParm)
%    fprintsf('Supplied ParmMask of nonmatching size');
%    ParmMask = false(size(InitParm));
%end

%if nargin < 5
%    placement = 'right';
%end

%fixedparams=[0 20];
switchxy=isa(Data,'dip_image');  % Dip images need coordinate switching of X and Y

if prod(ROISize) > 0
    MultiGaussSigmaMSE(double(Data),mymethod,numdims,double(centerzero),double(switchxy),[0 InitParm]',ROISize)  % submits the data, which will be copied inside the routine for later comparison with the simulation
else
    MultiGaussSigmaMSE(double(Data),mymethod,numdims,double(centerzero),double(switchxy))  % submits the data, which will be copied inside the routine for later comparison with the simulation
end
%[res, fitted] = MultiGaussSigmaMSE(InitParm');
%params=fixedparams;

tic
if (0)  % old method
    options=optimset('Display','notify','TolX',10^-4,'TolFun',10^-8,'MaxIter',maxiter); % ,'Method','newton0');
    %[params,msevalue]=fminsearch(@MultiGaussSigmaMSE,InitParm',options);
    params = fminsearch(@MultiGaussSigmaMSE,InitParm',options);
else    % new method with smarter amd faster optimisation routine from http://www.cs.ubc.ca/~schmidtm
    %options = struct('Display','on','notify',1,'numDiff',1,'TolX',10^-4,'TolFun',10^-8,'MaxIter',maxiter);
    options = struct('Display','on','notify',1,'TolX',10^-5,'TolFun',10^-7,'MaxIter',maxiter,'Method','pcg');
    %[params,msevalue,moreinfo]=minFunc(@(P)MultiGaussSigmaMSE(ParmMask' .* InitParm' + ~ParmMask' .* P),InitParm',options);
    [params,msevalue,moreinfo]=minFunc(@myMultiGaussSigmaMSE,tovary,options);
    %params = minFunc(@myMultiGaussSigmaMSE, tovary, options);
end
toc
params = ReconstructParam(params);
[res, mygrad, fitted,residual] = MultiGaussSigmaMSE(params);  % to generate the final result image and residual
params = params';
fitted=dip_image(permute(fitted,[2 1 3]));
residual=dip_image(permute(residual,[2 1 3]));
end

function [res,mygrad] = myMultiGaussSigmaMSE(tovary)
    param = ReconstructParam(tovary);
    if nargout>1
        [res,mygrad] = MultiGaussSigmaMSE(param);
        mygrad = CondenseParam(mygrad);
    else
        res = MultiGaussSigmaMSE(param);
    end
    if (1)
        % param
        [res2, mygrad2, fitted,residual] = MultiGaussSigmaMSE(param);  % to generate the final result image and residual
        dipshow(2,dip_image(residual));
        dipshow(3,dip_image(fitted));
        drawnow;
    end
end

% The function below condenses a set of full parameters and extracts only the ones to vary
function tovary = CondenseParam(FullParams)
    global ParmVarIdxFull
    global GlobVarBack
    global ParmGlobIdx   % not ParmGlobIdxFull, which contains all the global entries for expansion, but only the first ones for contraction
    if GlobVarBack == 1
        tovary=FullParams([1;1+ParmGlobIdx';1+ParmVarIdxFull]);
    else
        tovary=FullParams([1+ParmGlobIdx';1+ParmVarIdxFull]);
    end
end

% The function below takes a set of variable parameters and reconstructs the full parameter set
function param = ReconstructParam(InParam)
    global ParmVarIdxFull
    global ParmConstIdxFull
    global ParmGlobIdxFull
    global ParmConsts
    global nrspots
    global GlobVarBack

    param = zeros(1 + length(ParmVarIdxFull) + length(ParmConstIdxFull) + length(ParmGlobIdxFull), 1);
    if GlobVarBack == 1
        param(1) = InParam(1);
        InParam(1) = [];
    else
        param(1) = ParmConsts(1);
    end
    gparam = InParam(1 : (length(ParmGlobIdxFull) / nrspots));
    if ~isempty(gparam)
        InParam(1 : length(gparam)) = [];
        gparam = repmat(gparam, nrspots, 1);
    end
    gparam = gparam(:);
    param(ParmGlobIdxFull + 1) = gparam;
    param(ParmVarIdxFull + 1) = InParam;
    param(ParmConstIdxFull + 1) = ParmConsts((2 - GlobVarBack) : end);
end