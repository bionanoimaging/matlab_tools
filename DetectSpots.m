% [SpotList,Bg]=DetectSpots(InputImg,thresh,threshres,startsigma,method,ROISize,MaxSpotSize,MinInt)   : detects spots
% thresh : threshold in image to place new spots
% threshres : threshold in residuum to place new spots
% startsigma : starting sigma for fitting
% method : e.g 'idiv' method to use while fitting, see FitDataNDFastMask
% ROISize : X and Y radius of the ROI 
% MaxSpotSize : maximal size allowed for spots
% MinInt : minimu intensity allowed for spots


function [SpotList,Bg]=DetectSpots(InputImg,thresh,threshres,startsigma,method,ROISize,MaxSpotSize,MinInt)

maxiter=500;

if nargin < 5
    method='idiv';
end

if nargin < 6
    ROISize=[18 18];
end

if nargin < 7
    MaxSpotSize=20;
end

if nargin < 8
    MinInt=10;
end

filtered=gaussf(InputImg,max(6,startsigma));
Bg = min(filtered);

filtered=gaussf(InputImg,max(4,startsigma));
filteredorig=filtered;
SpotList=[];
n=1;
BgFlag=1;   % start with a fixed Background
while 1
    spots=dip_maxima(filtered,[],1,1);
    spots=((spots .* filtered) > thresh);  % In first round image. Later the residuals
    
    meas=measure(spots,filteredorig,{'Mean','Gravity','MaxVal'});
    if isempty(meas)
        fprintf('Found no more spots\n');
        break;
    end
    
    NewSpots = [meas.MaxVal-Bg; meas.Gravity(1,:);meas.Gravity(2,:); repmat(startsigma,[1 length(meas)]);repmat(startsigma,[1 length(meas)])];
    numSpots=size(NewSpots,2);
    fprintf('Step %d, adding %d spots\n',n,numSpots);
    tic
    initParms = [Bg,SpotList,NewSpots(:)'];
    % initParms=initParms(1:11)
    [params, res, fitted, residual] = FitDataNDFastMask(initParms, [BgFlag 0 0 0 0 0],InputImg, 2, maxiter, method, 0, ROISize);
    toc
    Bg=params(1);
    SpotList=params(2:end);
    
    SpotList=reshape(SpotList,[5,length(SpotList)/5]);
    numSpots=size(SpotList,2);
    SpotList(:,find(SpotList(5,:) > MaxSpotSize))=[];  % Delete wrong spots
    SpotList(:,find(SpotList(4,:) > MaxSpotSize))=[];  % Delete wrong spots
    SpotList(:,find(SpotList(5,:) < 0))=[];  % Delete wrong spots
    SpotList(:,find(SpotList(4,:) < 0))=[];  % Delete wrong spots
    SpotList(:,find(SpotList(1,:) < MinInt))=[];  % Delete wrong spots
    fprintf('Deleted %d spots\n',numSpots-size(SpotList,2))
    SpotList=SpotList(:)';
    
    initParms= [Bg,SpotList];
    % BgFlag=0;  % from the second run on, the background is fitted as well
    [params, res, fitted2, residual] = FitDataNDFastMask(initParms, [BgFlag 0 0 0 0 0],InputImg, 2, maxiter, method, 0, ROISize);
    Bg=params(1);
    SpotList=params(2:end);

    residual=InputImg-fitted2; % To make it work for ROIs as well
    % filtered = gaussf(abs(residual),4);
    filtered = gaussf(residual,max(2,startsigma));
    n=n+1;
    thresh=threshres;
end
SpotList=reshape(params(2:end),5,(size(params,2)-1)/5);

fprintf('Number of spots found: %d\n',size(SpotList,2))

% tmp=SpotList(2,:);   % Swap X and Y coordinates of position
% SpotList(2,:)=SpotList(3,:);
% SpotList(3,:)=tmp;
% 
% tmp=SpotList(4,:);   % Swap X and Y coordinates of sizes
% SpotList(4,:)=SpotList(5,:);
% SpotList(5,:)=tmp;
