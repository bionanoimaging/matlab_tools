% FitDataNDFast : Fits several ND Gaussians to N-dimensional image data
% Synopsis: [fitted,params,myfunct] = FitDataNDFast(InitParm,Data,numdims,maxiter,method, fixedparams)
%
% InitParm : Vector of initial parameters first row is global parameters,
% meaning is as follows: [global background, global slope X, global slope Y, .., posy1, posx1, .., width X1, width Y1, .., intensity1, posx2, posy2, ...]
% other rows are local parameters
% Data : Experimental multidimensional (e.g. image) data to be fitted
% maxiter : maximum number of iterations in fit (default = 300)
% method : figure of merit to use for fitting 'mes', 'idiv' or 'fidiv'
% fixedparams : Allows the user to fix any parameters to be either constant or globally fitted
%   it is a vector of the full length of standard parametervector with entries as follows:
%   0: fixed value. The value of the initialisation will be retained
%   1: variable value. This parameter will be fitted individually
%   2,3,..  : All corresponding parameters will be bound together with only one variable
%  default: []   (empty vector means all parameters are fitted)
%
%  Function is:   offset + sx*x + sy*y + brightness * exp(-( (x-x0)/wx)^2 + ...)
% FWHM is obtained by 2*sqrt(log(2))*wx
%
% Example:
% a=13+51.2*exp(-(((xx(20,20,10)-3.2)/2.2).^2+((yy(20,20,10)+2.8)/1.7).^2+(zz(20,20,10)/2.6).^2))
% an=noise(a,'poisson')
% [params,res,fitted,residual]=FitDataNDFast([10 0 0 0 1 1 1 3 3 2 40],an,3,300,'mse')
% fitted is the residual image (depending on method)
% params contains the result of the fit

% [params,res,fitted,residual]=FitDataNDFast([30 15 40 2 2],a,2,300,'idiv')

function [params,res,fitted,residual] = FitDataNDFast(InitParm,Data,numdims,maxiter,mymethod,fixedparams)

if nargin < 3
    numdims = 2;
end
if nargin < 4
    maxiter = 3000;
end

if nargin < 5
    mymethod = 'mse';
end
if nargin < 6
    fixedparams = [];
else
    if ~isempty(fixedparams)
        if norm(size(fixedparams)-size(InitParm)) > 0
            error('FitDataNDFast: Sizes of Initialisation and fixedparams have to be identical');
        end
    end
end

%if nargin < 5
%    placement = 'right';
%end

%fixedparams=[0 20];
%MultiGaussMSE(double(Data),mymethod,numdims)  % upload the data
MultiGaussMSEDeriv(double(Data),mymethod,numdims)  % upload the data, but account for the different meaning of X and Y in Matlab

InitParm=GaussArgPermuteXY(InitParm,numdims);
%[res, fitted] = MultiGaussMSE(InitParm');
[res, deriv, fitted] = MultiGaussMSEDeriv(InitParm');
%params=fixedparams;


if (0)  % old method
options=optimset('Display','notify','TolX',10^-4,'TolFun',10^-8,'MaxIter',maxiter);
[params,msevalue,fwd,resid]=fminsearch(@MultiGaussMSEDeriv,InitParm',options);
elseif (0)    % new method with smarter amd faster optimisation routine from http://www.cs.ubc.ca/~schmidtm
    options=struct('Display','off','notify',1,'numDiff',1,'TolX',10^-4,'TolFun',10^-8,'MaxIter',maxiter);
    [params,f,exitflag,output]=minFunc(@MultiGaussMSEDeriv,InitParm',options);
else
    options=struct('Display','off','notify',1,'TolX',10^-4,'TolFun',10^-8,'MaxIter',maxiter);
    if isempty(fixedparams)
        [params,f,exitflag,output]=minFunc(@MultiGaussMSEDeriv,InitParm',options);
    else
        [InitParm,constvec] = CondenseInitVecWithFixedVec(InitParm,fixedparams);
        f = @(fixedvec,constvec,numdims,x)( @(x)(FixedMultiGaussDeriv(x,fixedvec,constvec,numdims)));  % binds the other variables to constants
        [params,f,exitflag,output]=minFunc(f(fixedparams,constvec,numdims),InitParm',options);
        params=PackWithFixedVec(params,fixedparams,constvec,numdims);
    end
end

%[res, fitted,residual] = MultiGaussMSE(params);
[res, deriv, fitted,residual] = MultiGaussMSEDeriv(params);
params = GaussArgPermuteXY(params',numdims);
fitted=mat2im(fitted);
residual=mat2im(residual);
